#include "lmspredicter.h"

#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

LMSPredicter::~LMSPredicter()
{
	if(_readThread != 0)
		_readThread->join();
	clearBuffers();
}

void LMSPredicter::clearBuffers()
{
	for(MC2x2* buffer : _buffers)
		delete[] buffer;
}

void LMSPredicter::Start()
{
	boost::mutex::scoped_lock lock(_mutex);
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ROArrayColumn<casa::Complex> dataColumn(_ms, _ms.columnName(casa::MSMainEnums::DATA));
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	_bandData = BandData(_ms.spectralWindow());
	
	if(!_useModelColumn)
	{
		casa::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
		casa::MEpoch startTime = timeColumn(_startRow);
		size_t antennaCount = _ms.antenna().nrow();
		_dftInput->InitializeBeamBuffers(antennaCount, _bandData.ChannelCount());
		_predicter.reset(new DFTPredictionAlgorithm(*_dftInput, _bandData));
		if(_applyBeam)
		{
			_beamEvaluator.reset(new LBeamEvaluator(_ms));
			_beamEvaluator->SetTime(startTime);
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
		_workThreadGroup.reset(new boost::thread_group());
	_readThread.reset(new boost::thread(&LMSPredicter::ReadThreadFunc, this));
}
	
void LMSPredicter::ReadThreadFunc()
{
	size_t actualThreadCount = _threadCount;
	if(!_useModelColumn)
	{
		if(_dftInput->ComponentCount() == 0)
			actualThreadCount = 1;
		for(size_t i=0; i!=actualThreadCount; ++i)
			_workThreadGroup->add_thread(new boost::thread(&LMSPredicter::PredictThreadFunc, this));
	}
	
	boost::mutex::scoped_lock lock(_mutex);

	casa::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casa::MSMainEnums::UVW));
	casa::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
	std::unique_ptr<casa::ROArrayColumn<casa::Complex>> modelColumn;
	
	RowData rowData;
	
	casa::Array<casa::Complex> modelData;
	if(_useModelColumn)
	{
		modelColumn.reset(new casa::ROArrayColumn<casa::Complex>(_ms, _ms.columnName(casa::MSMainEnums::MODEL_DATA)));
		modelData = casa::Array<casa::Complex>(modelColumn->shape(0));
	}

	size_t timeIndex = 0;
	casa::MEpoch previousTime = timeColumn(_startRow);
	for(size_t rowIndex=_startRow; rowIndex!=_endRow; ++rowIndex)
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
				// Stop all threads, then update beam values, then restart threads.
				_workLane.write_end();
				_workThreadGroup->join_all();
				
				_workLane.clear();
				_beamEvaluator->SetTime(time);
				_predicter->UpdateBeam(*_beamEvaluator);
				_workThreadGroup.reset(new boost::thread_group());
				for(size_t i=0; i!=actualThreadCount; ++i)
					_workThreadGroup->add_thread(new boost::thread(&LMSPredicter::PredictThreadFunc, this));
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
				casa::Complex* inptr = modelData.cbegin();
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
		_workThreadGroup->join_all();
	}
	_outputLane.write_end();
}

void LMSPredicter::PredictThreadFunc()
{
	RowData rowData;
	while(_workLane.read(rowData))
	{
		MC2x2 *valIter = rowData.modelData;
		for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
		{
			double lambda = _bandData.ChannelWavelength(ch);
			_predicter->Predict(*valIter, rowData.u/lambda, rowData.v/lambda, rowData.w/lambda, ch, rowData.a1, rowData.a2);
			++valIter;
		}
		_outputLane.write(rowData);
	}
}
