#include "lmspredicter.h"

#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <boost/bind.hpp>

LMSPredicter::~LMSPredicter()
{
	if(_readThread != nullptr)
		_readThread->join();
	clearBuffers();
}

void LMSPredicter::InitializeInput(const Model& model)
{
	casacore::MSField fTable(_ms.field());
	casacore::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	casacore::MDirection refDir = refDirColumn(0);
	casacore::Vector<casacore::Double> val = refDir.getValue().get();
	double ra = val[0];
	double dec = val[1];
	
	BandData band(_ms.spectralWindow());
	
	_dftInput.InitializeFromModel(model, ra, dec, band);
}

void LMSPredicter::clearBuffers()
{
	for(MC2x2* buffer : _buffers)
		delete[] buffer;
}

void LMSPredicter::Start()
{
	std::lock_guard<std::mutex> lock(_mutex);
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casacore::ROArrayColumn<casacore::Complex> dataColumn(_ms, _ms.columnName(casacore::MSMainEnums::DATA));
	casacore::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	_bandData = BandData(_ms.spectralWindow());
	if(_endChannel == 0)
	{
		_startChannel = 0;
		_endChannel = _bandData.ChannelCount();
	}
	
	if(!_useModelColumn)
	{
		casacore::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
		casacore::MEpoch startTime = timeColumn(_startRow);
		size_t antennaCount = _ms.antenna().nrow();
		_dftInput.InitializeBeamBuffers(antennaCount, _bandData.ChannelCount());
		_predicter.reset(new DFTPredictionAlgorithm(_dftInput, _bandData));
		if(_applyBeam)
		{
			_beamEvaluator.reset(new LBeamEvaluator(_ms));
			_beamEvaluator->SetTime(startTime);
			_predicter->UpdateBeam(*_beamEvaluator, _startChannel, _endChannel);
		}
		else {
			_dftInput.SetUnitaryBeam();
		}
	}
	
	// Create buffers
	if(_buffers.empty())
	{
		_buffers.resize(_laneSize);
		for(size_t i=0; i!=_laneSize; ++i)
		{
			RowData rowData;
			_buffers[i] = new MC2x2[_bandData.ChannelCount()];
			rowData.modelData = _buffers[i];
			_availableBufferLane.write(rowData);
		}
	}

	// Start all threads
	if(!_useModelColumn)
		_workThreadGroup.clear();
	_readThread.reset(new std::thread(&LMSPredicter::ReadThreadFunc, this));
}
	
void LMSPredicter::ReadThreadFunc()
{
	size_t actualThreadCount = _threadCount;
	boost::asio::io_service::work work(_ioService);
	if(!_useModelColumn)
	{
		if(_dftInput->ComponentCount() == 0)
			actualThreadCount = 1;
		for(size_t i=0; i!=actualThreadCount; ++i)
		{
			_workThreadGroup.emplace_back(boost::bind(&boost::asio::io_service::run, &_ioService));
			_ioService.post(boost::bind(&LMSPredicter::PredictThreadFunc, this));
		}
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
			// If we get to a next timestep and are applying beam, update beam.
			if(!_useModelColumn && _applyBeam && _beamEvaluator->Time().getValue() != time.getValue())
			{
				// Put all threads on wait, then update beam values, then restart threads.
				_workLane.write_end();
				_barrier.wait();
				
				_workLane.clear();
				_beamEvaluator->SetTime(time);
				_predicter->UpdateBeam(*_beamEvaluator, _startChannel, _endChannel);
				for(size_t i=0; i!=actualThreadCount; ++i)
					_ioService.post(boost::bind(&LMSPredicter::PredictThreadFunc, this));
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
				MC2x2 *outptr = rowData.modelData;
				casacore::Complex* inptr = modelData.cbegin();
				for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
				{
					for(size_t p=0; p!=4; ++p)
					{
						outptr->Data()[p] = *inptr;
						++inptr;
					}
					++outptr;
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
		_barrier.wait();
		_ioService.stop();
		for(std::thread& t : _workThreadGroup)
			t.join();
		_workThreadGroup.clear();
	}
	_outputLane.write_end();
}

void LMSPredicter::PredictThreadFunc()
{
	RowData rowData;
	while(_workLane.read(rowData))
	{
		MC2x2 *valIter = rowData.modelData;
		for(size_t ch=_startChannel; ch!=_endChannel; ++ch)
		{
			double lambda = _bandData.ChannelWavelength(ch);
			_predicter->Predict(*valIter, rowData.u/lambda, rowData.v/lambda, rowData.w/lambda, ch, rowData.a1, rowData.a2);
			++valIter;
		}
		_outputLane.write(rowData);
	}
	_barrier.wait();
}
