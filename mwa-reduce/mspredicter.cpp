#include "mspredicter.h"

MSPredicter::~MSPredicter()
{
	_readThread->join();
	clearBuffers();
}

void MSPredicter::clearBuffers()
{
	for(std::vector<std::complex<double>*>::iterator i=_buffers.begin(); i!=_buffers.end(); ++i)
		delete[] *i;
}

void MSPredicter::Start()
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
	_predicter->Initialize(_model, &_beamEvaluator);
	
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
	size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	if(_model.SourceCount() == 0)
		cpuCount = 1;
	for(size_t i=0; i!=cpuCount; ++i)
		_workThreadGroup.add_thread(new boost::thread(&MSPredicter::PredictThreadFunc, this));
	_readThread.reset(new boost::thread(&MSPredicter::ReadThreadFunc, this));
}
	
void MSPredicter::ReadThreadFunc()
{
	boost::mutex::scoped_lock lock(_mutex);

	casa::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casa::MSMainEnums::UVW));
	
	RowData rowData;

	for(size_t rowIndex=0; rowIndex!=_ms.nrow(); ++rowIndex)
	{
		if(ant1Column(rowIndex) != ant2Column(rowIndex))
		{
			casa::Array<double> uvwArray = uvwColumn(rowIndex);
			casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
			double u = *uvwI; ++uvwI;
			double v = *uvwI; ++uvwI;
			double w = *uvwI;
			lock.unlock();
			
			_availableBufferLane.read(rowData);
			rowData.u = u;
			rowData.v = v;
			rowData.w = w;
			rowData.rowIndex = rowIndex;
			_workLane.write(rowData);
			
			lock.lock();
		}
	}
	
	lock.unlock();
	_workLane.write_end();
	_workThreadGroup.join_all();
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
			_predicter->Predict4(valIter, _model, rowData.u/lambda, rowData.v/lambda, rowData.w/lambda, ch);
			valIter += 4;
		}
		_outputLane.write(rowData);
	}
}
