#include "mspredicter.h"

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

MSPredicter::~MSPredicter()
{
	if(_readThread != 0)
		_readThread->join();
	clearBuffers();
}

void MSPredicter::clearBuffers()
{
	for(std::vector<std::complex<double>*>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		delete[] *i;
}

void MSPredicter::Start(bool reportSources)
{
	std::lock_guard<std::mutex> lock(_mutex);
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casacore::ROArrayColumn<casacore::Complex> dataColumn(_ms, _ms.columnName(casacore::MSMainEnums::DATA));
	casacore::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	_bandData.reset(new BandData(_ms.spectralWindow()));
	_channelCount = _bandData->ChannelCount();
	
	casacore::MSField fieldTable = _ms.field();
	casacore::ROArrayColumn<double> phaseDirColumn(fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	if(phaseDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casacore::Array<double> phaseDir = phaseDirColumn(0);
	casacore::Array<double>::const_iterator phaseDirIter = phaseDir.begin();
	long double phaseCentreRA = *phaseDirIter; ++phaseDirIter;
	long double phaseCentreDec = *phaseDirIter;
	// By setting the time beforehand, we don't waste time calculating a time step we don't need.
	if(!_useModelColumn)
	{
		casacore::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
		casacore::MEpoch startTime = timeColumn(_startRow);
		_predicter.reset(new Predicter(phaseCentreRA, phaseCentreDec, _bandData->LowestFrequency(), _bandData->HighestFrequency(), _channelCount));
		if(_applyBeam)
		{
			_beamEvaluator.reset(new BeamEvaluator(_ms, false, _mwaPath));
			_beamEvaluator->SetTime(startTime);
			_predicter->Initialize(_model, _solutionFile, &*_beamEvaluator);
		}
		else
			_predicter->Initialize(_model, _solutionFile);
		if(reportSources)
			_predicter->ReportSources(_model);
	}
	
	// Create buffers
	if(_buffers.empty())
	{
		_buffers.resize(_laneSize);
		for(size_t i=0; i!=_laneSize; ++i)
		{
			RowData rowData;
			_buffers[i] = new std::complex<double>[_channelCount*4];
			rowData.modelData = _buffers[i];
			_availableBufferLane.write(rowData);
		}
	}

	// Start all threads
	if(!_useModelColumn)
		_workThreadGroup.clear();
	_readThread.reset(new std::thread(&MSPredicter::ReadThreadFunc, this));
}
	
void MSPredicter::ReadThreadFunc()
{
	size_t actualThreadCount = _threadCount;
	if(!_useModelColumn)
	{
		if(_model.SourceCount() == 0)
			actualThreadCount = 1;
		for(size_t i=0; i!=actualThreadCount; ++i)
			_workThreadGroup.emplace_back(&MSPredicter::PredictThreadFunc, this);
	}
	
	std::unique_lock<std::mutex> lock(_mutex);

	casacore::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casacore::MSMainEnums::UVW));
	casacore::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
	std::unique_ptr<casacore::ROArrayColumn<casacore::Complex>> modelColumn;
	
	RowData rowData;
	
	casacore::Array<casacore::Complex> modelData;
	if(_useModelColumn)
	{
		modelColumn.reset(new casacore::ROArrayColumn<casacore::Complex>(_ms, _ms.columnName(casacore::MSMainEnums::MODEL_DATA)));
		modelData = casacore::Array<casacore::Complex>(modelColumn->shape(0));
	}

	size_t timeIndex = 0;
	casacore::MEpoch previousTime = timeColumn(_startRow);
	for(size_t rowIndex=_startRow; rowIndex!=_endRow; ++rowIndex)
	{
		size_t
			a1 = ant1Column(rowIndex),
			a2 = ant2Column(rowIndex);
		casacore::MEpoch time = timeColumn(rowIndex);
		if(a1 != a2)
		{
			casacore::Array<double> uvwArray = uvwColumn(rowIndex);
			casacore::Array<double>::const_contiter uvwI = uvwArray.cbegin();
			double u = *uvwI; ++uvwI;
			double v = *uvwI; ++uvwI;
			double w = *uvwI;
			if(_useModelColumn)
				modelColumn->get(rowIndex, modelData);
			lock.unlock();
			
			if(time.getValue() != previousTime.getValue())
			{
				++timeIndex;
				previousTime = time;
			}
			if(!_useModelColumn && _applyBeam && _beamEvaluator->Time().getValue() != time.getValue())
			{
				// Stop all threads, then update beam values, then restart threads.
				_workLane.write_end();
				for(std::thread& t : _workThreadGroup)
					t.join();
				
				_workLane.clear();
				_beamEvaluator->SetTime(time);
				_predicter->UpdateBeam(_model);
				_workThreadGroup.clear();
				for(size_t i=0; i!=actualThreadCount; ++i)
					_workThreadGroup.emplace_back(&MSPredicter::PredictThreadFunc, this);
			}
			
			_availableBufferLane.read(rowData);
			rowData.u = u;
			rowData.v = v;
			rowData.w = w;
			rowData.rowIndex = rowIndex;
			rowData.a1 = a1;
			rowData.a2 = a2;
			rowData.timeIndex = timeIndex;
			if(_useModelColumn)
			{
				std::complex<double> *outptr = rowData.modelData;
				casacore::Complex* inptr = modelData.cbegin();
				for(size_t ch=0; ch!=_channelCount*4; ++ch)
				{
					*outptr = *inptr;
					++outptr; ++inptr;
				}
				_outputLane.write(rowData);
			}
			else {
				_workLane.write(rowData);
			}
			
			lock.lock();
		}
	}
	
	lock.unlock();
	if(!_useModelColumn)
	{
		_workLane.write_end();
		for(std::thread& t : _workThreadGroup)
			t.join();
	}
	_outputLane.write_end();
}

void MSPredicter::PredictThreadFunc()
{
	RowData rowData;
	while(_workLane.read(rowData))
	{
		std::complex<double> *valIter = rowData.modelData;
		for(size_t ch=0; ch!=_channelCount; ++ch)
		{
			double lambda = _bandData->ChannelWavelength(ch);
			_predicter->Predict4(valIter, _model, rowData.u/lambda, rowData.v/lambda, rowData.w/lambda, ch, rowData.a1, rowData.a2);
			valIter += 4;
		}
		_outputLane.write(rowData);
	}
}
