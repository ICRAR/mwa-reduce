#include "mspredicter.h"

#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

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
	boost::mutex::scoped_lock lock(_mutex);
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ROArrayColumn<casa::Complex> dataColumn(_ms, _ms.columnName(casa::MSMainEnums::DATA));
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	_bandData.reset(new BandData(_ms.spectralWindow()));
	_channelCount = _bandData->ChannelCount();
	
	casa::MSField fieldTable = _ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	
	_predicter.reset(new Predicter(phaseCentreRA, phaseCentreDec, _bandData->LowestFrequency(), _bandData->HighestFrequency(), _channelCount));
	if(_applyBeam)
		_predicter->Initialize(_model, _solutionFile, &_beamEvaluator);
	else
		_predicter->Initialize(_model, _solutionFile);
	if(reportSources)
		_predicter->ReportSources(_model);
	
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
	_workThreadGroup.reset(new boost::thread_group());
	_readThread.reset(new boost::thread(&MSPredicter::ReadThreadFunc, this));
}
	
void MSPredicter::ReadThreadFunc()
{
	size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	if(_model.SourceCount() == 0)
		cpuCount = 1;
	for(size_t i=0; i!=cpuCount; ++i)
		_workThreadGroup->add_thread(new boost::thread(&MSPredicter::PredictThreadFunc, this));
	
	boost::mutex::scoped_lock lock(_mutex);

	casa::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casa::MSMainEnums::UVW));
	casa::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
	
	RowData rowData;

	for(size_t rowIndex=0; rowIndex!=_ms.nrow(); ++rowIndex)
	{
		size_t
			a1 = ant1Column(rowIndex),
			a2 = ant2Column(rowIndex);
		casa::MEpoch time = timeColumn(rowIndex);
		if(a1 != a2)
		{
			casa::Array<double> uvwArray = uvwColumn(rowIndex);
			casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
			double u = *uvwI; ++uvwI;
			double v = *uvwI; ++uvwI;
			double w = *uvwI;
			lock.unlock();
			
			if(_beamEvaluator.Time().getValue() != time.getValue())
			{
				// Stop all threads, then update beam values, then restart threads.
				_workLane.write_end();
				_workThreadGroup->join_all();
				
				_workLane.clear();
				_beamEvaluator.SetTime(time);
				_predicter->UpdateBeam(_model);
				_workThreadGroup.reset(new boost::thread_group());
				for(size_t i=0; i!=cpuCount; ++i)
					_workThreadGroup->add_thread(new boost::thread(&MSPredicter::PredictThreadFunc, this));
			}
			
			_availableBufferLane.read(rowData);
			rowData.u = u;
			rowData.v = v;
			rowData.w = w;
			rowData.rowIndex = rowIndex;
			rowData.a1 = a1;
			rowData.a2 = a2;
			_workLane.write(rowData);
			
			lock.lock();
		}
	}
	
	lock.unlock();
	_workLane.write_end();
	std::cout << "Joining all.\n";
	_workThreadGroup->join_all();
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
