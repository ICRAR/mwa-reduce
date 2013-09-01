#include "spectrumsubtractor.h"

#include "model.h"
#include "predicter.h"
#include "matrix2x2.h"

#include <ms/MeasurementSets/MeasurementSet.h>

const size_t SpectrumSubtractor::BUFFER_COUNT = 64;

SpectrumSubtractor::SpectrumSubtractor(casa::MeasurementSet& ms, Model& model) :
	_ms(ms),
	_model(model),
	_bandData(ms.spectralWindow()),
	_subtractWorkLane(BUFFER_COUNT),
	_subtractAvailableBufferLane(BUFFER_COUNT),
	_timestepCount(0),
	_fittingInterval(1)
{
	_cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
}

SpectrumSubtractor::~SpectrumSubtractor()
{
}

void SpectrumSubtractor::Perform()
{
	_dataColumn.reset(new casa::ArrayColumn<casa::Complex>(_ms, _dataColumnName));
	_antenna1Column.reset(new casa::ROScalarColumn<int> (_ms, _ms.columnName(casa::MSMainEnums::ANTENNA1)));
	_antenna2Column.reset(new casa::ROScalarColumn<int> (_ms, _ms.columnName(casa::MSMainEnums::ANTENNA2)));
	_uvwColumn.reset(new casa::ROArrayColumn<double> (_ms, _ms.columnName(casa::MSMainEnums::UVW)));
	casa::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(_ms, _ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<double> timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
	
	casa::IPosition dataShape = _dataColumn->shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Need 4 polarizations");
	
	double time = timeColumn(0);
	size_t intervalStartRow = 0;
	
	initSources();
	initPredictors();
	countTimesteps(timeColumn);
	size_t intervalCount = _timestepCount / _fittingInterval;
	if(intervalCount == 0) intervalCount = 1;
	std::cout << "Measuring spectra for " << intervalCount << " intervals...\n";
	
	_spectrumSums.assign(4 * _bandData.ChannelCount() * _model.SourceCount(), 0.0);
	_spectrumWeights.assign(4 * _bandData.ChannelCount() * _model.SourceCount(), 0.0);
	_totalFluxPerSource.assign(4 * _model.SourceCount(), 0.0);
	_totalFluxWeightPerSource.assign(4 * _model.SourceCount(), 0.0);
	
	_dataBuffers.resize(BUFFER_COUNT);
	_weightBuffers.resize(BUFFER_COUNT);
	_flagBuffers.resize(BUFFER_COUNT);
	for(size_t i=0; i!=BUFFER_COUNT; ++i)
	{
		_dataBuffers[i] = new casa::Array<casa::Complex>(dataShape);
		_weightBuffers[i] = new casa::Array<float>(dataShape);
		_flagBuffers[i] = new casa::Array<bool>(dataShape);
	}
	
	initMeasureThreadData();
	startMeasureThreads();
	
	size_t bufferIndex = 0, timeIndex = 0, intervalIndex = 0;
	for(size_t row=0; row!=_ms.nrow(); ++row)
	{
		size_t
			a1 = (*_antenna1Column)(row),
			a2 = (*_antenna2Column)(row);
		if(a1 != a2)
		{
			double thisTime = timeColumn(row);
			if(thisTime != time)
			{
				++timeIndex;
				time = thisTime;
				size_t nextIntervalTimestep = _timestepCount * (intervalIndex+1) / intervalCount;
				if(timeIndex == nextIntervalTimestep)
				{
					stopMeasureThreads();
					
					performSubtraction(intervalStartRow, row);
					
					intervalStartRow = row;
					startMeasureThreads();
					++intervalIndex;
				}
			}
			
			casa::Array<casa::Complex> &dataArray = *_dataBuffers[bufferIndex];
			casa::Array<float> &weightArray = *_weightBuffers[bufferIndex];
			casa::Array<bool> &flagArray = *_flagBuffers[bufferIndex];
				
			_dataColumn->get(row, dataArray);
			flagColumn.get(row, flagArray);
			weightColumn.get(row, weightArray);
			
			casa::Array<double> uvwArray = (*_uvwColumn)(row);
			casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
			double u = *uvwI; ++uvwI;
			double v = *uvwI; ++uvwI;
			double w = *uvwI;
				
			// Push data per source to threads
			for(size_t s=0; s!=_sources.size(); ++s)
			{
				size_t thread = s % _cpuCount;
				MeasureThreadInfo task;
				task.data = dataArray.cbegin();
				task.weights = weightArray.cbegin();
				task.flags = flagArray.cbegin();
				task.sourceIndex = s;
				task.u = u; task.v = v; task.w = w;
				task.a1 = a1; task.a2 = a2;
				_taskLanes[thread]->write(task);
			}
			
			bufferIndex = (bufferIndex + 1) % BUFFER_COUNT;
		}
	}
	stopMeasureThreads();
	performSubtraction(intervalStartRow, _ms.nrow());
	
	// Restore the model
	std::cout << "Subtracted absolute fluxes:\n";
	for(size_t s=0; s!=_model.SourceCount(); ++s)
	{
		ModelSource &source = _model.Source(s);
		double flux[4];
		std::cout << source.Name() << '\t';
		for(size_t p=0; p!=4; ++p)
		{
			flux[p] = _totalFluxPerSource[s*4 + p] / _totalFluxWeightPerSource[s*4 + p];
			if(p == 3)
				std::cout << flux[p] << '\n';
			else
				std::cout << flux[p] << '\t';
		}
		source.SetConstantTotalFlux(flux, _bandData.CentreFrequency());
	}
	
	_dataColumn.reset();
	_antenna1Column.reset();
	_antenna2Column.reset();
	_uvwColumn.reset();
}

void SpectrumSubtractor::initMeasureThreadData()
{
	_taskLanes.resize(_cpuCount);
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		// The task lane must issue a "wait" when its size is
		// BUFFER_COUNT-2, because the pushing thread can start
		// overwritting the n-1 th buffer while the popping thread
		// can just get the n-1 th buffer out.
		_taskLanes[i] = new lane<MeasureThreadInfo>(BUFFER_COUNT-2);
	}
}

void SpectrumSubtractor::initPredictors()
{
	casa::MSField fieldTable = _ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	
	for(std::vector<ModelSource>::iterator sourceIter=_sources.begin();
			sourceIter!=_sources.end(); ++sourceIter)
	{
		_predicters.push_back(std::unique_ptr<Predicter>(
			new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), _bandData.ChannelCount())));
		(*_predicters.rbegin())->Initialize(*sourceIter);
	}
}

void SpectrumSubtractor::countTimesteps(casa::ROScalarColumn<double>& timeColumn)
{
	double currentTimestep = timeColumn(0);
	_timestepCount = 1;
	for(size_t i=0; i!=_ms.nrow(); ++i)
	{
		if(timeColumn(i) != currentTimestep)
		{
			++_timestepCount;
			currentTimestep = timeColumn(i);
		}
	}
}

void SpectrumSubtractor::startMeasureThreads()
{
	_threadGroup.reset(new boost::thread_group());
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		_taskLanes[i]->clear();
		_threadGroup->add_thread(new boost::thread(&SpectrumSubtractor::measureThreadFunc, this, i));
	}
}

void SpectrumSubtractor::stopMeasureThreads()
{
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		_taskLanes[i]->write_end();
	}
	_threadGroup->join_all();
	_threadGroup.reset();
}

void SpectrumSubtractor::measureThreadFunc(size_t threadIndex)
{
	size_t channelCount = _bandData.ChannelCount();
	MeasureThreadInfo taskInfo;
	lane<MeasureThreadInfo> &taskLane = *_taskLanes[threadIndex];
	while(taskLane.read(taskInfo))
	{
		std::complex<float>* dataPtr = taskInfo.data;
		float* weightPtr = taskInfo.weights;
		bool* flagPtr = taskInfo.flags;
		
		size_t s = taskInfo.sourceIndex;
		double* measSumIter = &_spectrumSums[s * channelCount * 4];
		double* measWeightIter = &_spectrumWeights[s * channelCount * 4];
		
		Predicter& predicter = *_predicters[s];
		double
			u = taskInfo.u,
			v = taskInfo.v,
			w = taskInfo.w;
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double lambda = _bandData.ChannelWavelength(ch);
			Predicter::CNumType predicted[4];
			predicter.Predict4(predicted, _sources[s], u/lambda, v/lambda, w/lambda, ch, taskInfo.a1, taskInfo.a2);
			for(size_t p=0; p!=4; ++p)
			{
				double 
					real = dataPtr->real(),
					imag = dataPtr->imag();
				if(!(*flagPtr) && std::isfinite(real) && std::isfinite(imag))
				{
					// add real(data * conj(predicted))
					*measSumIter += *weightPtr * (real * predicted[0].real() + imag * predicted[0].imag());
					*measWeightIter += *weightPtr;
				}
				++dataPtr;
				++weightPtr;
				++flagPtr;
				++measSumIter;
				++measWeightIter;
			}
		}
	}
}

void SpectrumSubtractor::initSources()
{
	// Set all sources to flux 1 Jy
	for(Model::iterator sourceIter=_model.begin(); sourceIter!=_model.end(); ++sourceIter)
	{
		ModelSource newSource = *sourceIter;
		newSource.MakeUnitFlux();
		_sources.push_back(newSource);
	}
}

void SpectrumSubtractor::performSubtraction(size_t startRow, size_t endRow)
{
	// Normalize the spectrum
	double totalFlux = 0.0, totalWeight = 0.0;
	for(size_t s=0; s!=_sources.size(); ++s)
	{
		for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
		{
			size_t chIndex = (s*_bandData.ChannelCount() + ch) * 4;
			for(size_t p=0; p!=4; ++p)
			{
				size_t index = chIndex + p;
				double sum = _spectrumSums[index], weight = _spectrumWeights[index];
				if(std::isfinite(sum))
				{
					_totalFluxPerSource[s*4 + p] += sum;
					_totalFluxWeightPerSource[s*4 + p] += weight;
					
					totalFlux += sum;
					totalWeight += weight;
				}
				_spectrumSums[index] /= weight;
			}
		}
	}
	std::cout << "Flux: " << (2.0*totalFlux/totalWeight) << "\tWeight:" << totalWeight << '\n';
	
	_subtractWorkLane.clear();
	_subtractAvailableBufferLane.clear();
	
	for(size_t i=0; i!=BUFFER_COUNT; ++i)
	{
		SubtractThreadInfo info;
		info.readyForWrite = false;
		info.data = _dataBuffers[i];
		_subtractAvailableBufferLane.write(info);
	}
	
	startSubtractionThreads();
	
	for(size_t row=startRow; row!=endRow; ++row)
	{
		size_t a1 = (*_antenna1Column)(row), a2 = (*_antenna2Column)(row);
		if(a1 != a2)
		{
			SubtractThreadInfo info;
			_subtractAvailableBufferLane.read(info);
			if(info.readyForWrite)
			{
				_dataColumn->put(info.rowIndex, *info.data);
			}
			
			info.readyForWrite = false;
			info.rowIndex = row;
			_dataColumn->get(row, *info.data);
			
			casa::Array<double> uvwArray = (*_uvwColumn)(row);
			casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
			info.u = *uvwI; ++uvwI;
			info.v = *uvwI; ++uvwI;
			info.w = *uvwI;
			info.a1 = a1;
			info.a2 = a2;
			
			_subtractWorkLane.write(info);
		}
	}
	
	stopSubtractionThreads();
	
	while(!_subtractAvailableBufferLane.empty())
	{
		SubtractThreadInfo info;
		_subtractAvailableBufferLane.read(info);
		if(info.readyForWrite)
		{
			_dataColumn->put(info.rowIndex, *info.data);
		}
	}

	_spectrumSums.assign(4 * _bandData.ChannelCount() * _model.SourceCount(), 0.0);
	_spectrumWeights.assign(4 * _bandData.ChannelCount() * _model.SourceCount(), 0.0);
}

void SpectrumSubtractor::startSubtractionThreads()
{
	_threadGroup.reset(new boost::thread_group());
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		_threadGroup->add_thread(new boost::thread(&SpectrumSubtractor::subtractionThreadFunc, this));
	}
}

void SpectrumSubtractor::stopSubtractionThreads()
{
	_subtractWorkLane.write_end();
	_threadGroup->join_all();
	_threadGroup.reset();
}

void SpectrumSubtractor::subtractionThreadFunc()
{
	const size_t channelCount = _bandData.ChannelCount();
	SubtractThreadInfo info;
	while(_subtractWorkLane.read(info))
	{
		casa::Complex *data = info.data->cbegin();
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double lambda = _bandData.ChannelWavelength(ch);
			double u = info.u/lambda, v = info.v/lambda, w = info.w/lambda;
			std::complex<double> predictSum[4];
			for(size_t s=0; s!=_sources.size(); ++s)
			{
				Predicter &predicter = *_predicters[s];
				Predicter::CNumType predicted[4];
				predicter.Predict4(predicted, _sources[s], u, v, w, ch, info.a1, info.a2);
				for(size_t p=0; p!=4; ++p)
				{
					predicted[p] *= _spectrumSums[(s*channelCount + ch)*4 + p];
				}
				Matrix2x2::Add(predictSum, predicted);
			}
			for(size_t p=0; p!=4; ++p)
				data[ch*4 + p] -= casa::Complex(predictSum[p].real(), predictSum[p].imag());
		}
		
		info.readyForWrite = true;
		_subtractAvailableBufferLane.write(info);
	}
}
