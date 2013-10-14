#include "spectrummaker.h"

#include "mspredicter.h"
#include "banddata.h"
#include "beamevaluator.h"
#include "progressbar.h"
#include "modelsource.h"
#include "spectralenergydistribution.h"
#include "predicter.h"

#include <boost/thread/thread.hpp>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <limits>

SpectrumMaker::SpectrumMaker()
{
}

SpectrumMaker::~SpectrumMaker()
{
}

void SpectrumMaker::measure(const string& filename, const string& solutionsFile)
{
	casa::MeasurementSet ms(filename);
	
	/**
		* Read some meta data from the measurement set
		*/
	_bandData = BandData(ms.spectralWindow());
	size_t channelCount = _bandData.ChannelCount();
	
	casa::MSField fieldTable = ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	
	std::vector<std::complex<double>>
		measFlux(channelCount * _sources.size() * 4),
		measWeights(channelCount * _sources.size() * 4);
	
	MSPredicter modelPredicter(ms, _subtractedModel, solutionsFile);
	modelPredicter.SetApplyBeam(_applyBeam);
	
	std::vector<std::unique_ptr<Predicter>> predicters;
	for(std::vector<ModelSource>::iterator sourceIter=_sources.begin();
			sourceIter!=_sources.end(); ++sourceIter)
	{
		predicters.push_back(std::unique_ptr<Predicter>(
			new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), channelCount)));
		(*predicters.rbegin())->Initialize(*sourceIter);
	}
	
	_beamWeights[0].resize(_sources.size() * channelCount * 4);
	_beamWeights[1].resize(_sources.size() * channelCount * 4);
	if(_applyBeam)
		_beamEvaluator.reset(new BeamEvaluator(ms));
	
	const size_t BUFFER_COUNT = 16;
	size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	boost::thread_group threadGroup;
	std::vector<ao::lane<ThreadTaskInfo>*> taskLanes(cpuCount);
	for(size_t i=0; i!=cpuCount; ++i)
	{
		// The task lane must issue a "wait" when its size is
		// BUFFER_COUNT-2, because the pushing thread can start
		// overwritting the n-1 th buffer while the popping thread
		// can just get the n-1 th buffer out.
		taskLanes[i] = new ao::lane<ThreadTaskInfo>(BUFFER_COUNT-2);
		threadGroup.add_thread(new boost::thread(&SpectrumMaker::measureThreadFunc, this, &*taskLanes[i]));
	}
	
	std::vector<std::complex<double>*> dataBuffers(BUFFER_COUNT);
	std::vector<casa::Array<float>*> weightBuffers(BUFFER_COUNT);
	std::vector<casa::Array<bool>*> flagBuffers(BUFFER_COUNT);
	for(size_t i=0; i!=BUFFER_COUNT; ++i)
	{
		dataBuffers[i] = new std::complex<double>[channelCount * polarizationCount];
		weightBuffers[i] = new casa::Array<float>(dataShape);
		flagBuffers[i] = new casa::Array<bool>(dataShape);
	}
	casa::Array<casa::Complex> dataArray(dataShape);
	
	modelPredicter.Start();
	
	/**
		* Calculate spectra
		*/
	size_t bufferIndex = 0;
	MSPredicter::RowData rowData;
	size_t beamWeightIndex = 0;
	recalculateBeamWeights(beamWeightIndex);
	
	ProgressBar progress(std::string("Measure spectra in ") + filename);
	while(modelPredicter.GetNextRow(rowData))
	{
		size_t rowIndex = rowData.rowIndex;
		
		boost::mutex::scoped_lock lock(modelPredicter.IOMutex());
		progress.SetProgress(rowIndex, ms.nrow());
		
		// Cross correlation?
		if(rowData.a1 != rowData.a2)
		{
			std::complex<double> *data = dataBuffers[bufferIndex];
			casa::Array<float> &dataWeights = *weightBuffers[bufferIndex];
			casa::Array<bool> &flags = *flagBuffers[bufferIndex];
			
			dataColumn.get(rowIndex, dataArray);
			weightColumn.get(rowIndex, dataWeights);
			flagColumn.get(rowIndex, flags);
			casa::MEpoch time = timeColumn(rowIndex);
			lock.unlock();
			
			if(_applyBeam && time.getValue() != _beamEvaluator->Time().getValue())
			{
				std::cout << 'B' << std::flush;
				_beamEvaluator->SetTime(time);
				beamWeightIndex = (beamWeightIndex+1)%2;
				recalculateBeamWeights(beamWeightIndex);
			}
			
			casa::Array<casa::Complex>::const_contiter dataArrayIter = dataArray.cbegin();
			for(size_t i=0; i!=4*channelCount; ++i)
			{
				data[i] = std::complex<double>(dataArrayIter->real(), dataArrayIter->imag()) - rowData.modelData[i];
				++dataArrayIter;
			}
			
			for(size_t s=0; s!=_sources.size(); ++s)
			{
				size_t thread = s % cpuCount;
				ThreadTaskInfo task;
				task.flux = &measFlux[s * channelCount * 4];
				task.fluxWeights = &measWeights[s * channelCount * 4];
				task.data = data;
				task.dataWeights = dataWeights.cbegin();
				task.flags = flags.cbegin();
				task.u = rowData.u;
				task.v = rowData.v;
				task.w = rowData.w;
				task.beamWeights = &_beamWeights[beamWeightIndex][s * channelCount * 4];
				task.sourceIndex = s;
				task.predicter = &*predicters[s];
				task.a1 = rowData.a1;
				task.a2 = rowData.a2;
				taskLanes[thread]->write(task);
			}
			
			bufferIndex = (bufferIndex + 1) % BUFFER_COUNT;
		}
		
		modelPredicter.FinishRow(rowData);
	}
	
	for(size_t i=0; i!=cpuCount; ++i)
		taskLanes[i]->write_end();
	threadGroup.join_all();
	
	for(size_t i=0; i!=cpuCount; ++i)
		delete taskLanes[i];
	
	for(size_t i=0; i!=BUFFER_COUNT; ++i)
	{
		delete[] dataBuffers[i];
		delete weightBuffers[i];
		delete flagBuffers[i];
	}
	
	progress.SetProgress(ms.nrow(), ms.nrow());
	
	// Add to the total values
	for(size_t s=0; s!=_sources.size(); ++s)
	{
		Spectrum &spectrum = _spectrumPerSource[s];
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double freq = _bandData.ChannelFrequency(ch);
			spectrum.AddMeasurement(freq, &measFlux[(s * channelCount + ch) * 4], &measWeights[(s * channelCount + ch) * 4]);
		}
	}
}

void SpectrumMaker::recalculateBeamWeights(size_t beamWeightIndex)
{
	std::complex<double>* beamWeightPtr = &_beamWeights[beamWeightIndex][0];
	if(_applyBeam)
	{
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			const ModelSource& source = _sources[s];
			BeamEvaluator::PrecalcPosInfo posInfo;
			_beamEvaluator->PrecalculatePositionInfo(posInfo, source.Peak().PosRA(), source.Peak().PosDec());
			for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
			{
				_beamEvaluator->EvaluateAbsToApparentGain(posInfo, _bandData.ChannelFrequency(ch), beamWeightPtr);
				beamWeightPtr += 4;
			}
		}
	} else {
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
			{
				beamWeightPtr[0] = 1.0; beamWeightPtr[1] = 0.0;
				beamWeightPtr[2] = 0.0; beamWeightPtr[3] = 1.0;
				beamWeightPtr += 4;
			}
		}
	}
}

void SpectrumMaker::recalculateBeamWeightsThreadFunc(ao::lane<BeamEvalTaskInfo> *taskLane)
{
	BeamEvalTaskInfo info;
	while(taskLane->read(info))
	{
		std::complex<double>* weights = info.weights;
		const ModelSource& src = *info.source;
		for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
		{
			_beamEvaluator->EvaluateAbsToApparentGain(src.Peak().PosRA(), src.Peak().PosDec(), _bandData.ChannelFrequency(ch), weights);
			weights += 4;
		}
	}
}

void SpectrumMaker::measureThreadFunc(ao::lane<SpectrumMaker::ThreadTaskInfo>* taskLane)
{
	size_t channelCount = _bandData.ChannelCount();
	ThreadTaskInfo taskInfo;
	while(taskLane->read(taskInfo))
	{
		std::complex<double>* dataPtr = taskInfo.data;
		float* weightPtr = taskInfo.dataWeights;
		std::complex<double>* beamWeightPtr = taskInfo.beamWeights;
		bool* flagPtr = taskInfo.flags;
		std::complex<double>* measFluxIter = taskInfo.flux;
		std::complex<double>* measWeightIter = taskInfo.fluxWeights;
		Predicter& predicter = *taskInfo.predicter;
		double
			u = taskInfo.u,
			v = taskInfo.v,
			w = taskInfo.w;
		size_t s = taskInfo.sourceIndex;
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double lambda = _bandData.ChannelWavelength(ch);
			Predicter::CNumType predicted[4];
			predicter.Predict4(predicted, _sources[s], u/lambda, v/lambda, w/lambda, ch, taskInfo.a1, taskInfo.a2);
			bool sampleGood = true;
			double weightScalar = 0.0;
			for(size_t p=0; p!=4; ++p)
			{
				double 
					real = dataPtr->real(),
					imag = dataPtr->imag();
				if(*flagPtr || !std::isfinite(real) || !std::isfinite(imag))
					sampleGood = false;
				else
					weightScalar += *weightPtr;
				
				++weightPtr;
				++flagPtr;
			}
			
			std::complex<double> visSample[4];
			if(sampleGood) {
				Matrix2x2::ATimesHermB(visSample, dataPtr, predicted);
				Matrix2x2::PlusHermATimesB(visSample, dataPtr, predicted);
			}
			
			dataPtr += 4;
			
			if(sampleGood)
			{
				weightScalar = 0.25 * weightScalar;
			
				// Calculate Flux += w B* V B  (from: w (B* B) B^-1 V B*^-1 (B* B))
				// w = data weight, B = beam weight, V = vis
				std::complex<double> temp[4];
				Matrix2x2::HermATimesB(temp, beamWeightPtr, visSample);
				Matrix2x2::ScalarMultiply(temp, weightScalar * 0.5); // Divide factor of 2 because we add both normal and conjugate at once
				Matrix2x2::PlusATimesB(measFluxIter, temp, beamWeightPtr);
				//Matrix2x2::MultiplyAdd(measFluxIter, visSample, weightScalar);
				
				// Calculate Weight += w B* B
				Matrix2x2::HermATimesB(temp, beamWeightPtr, beamWeightPtr);
				std::complex<double> temp2[4];
				Matrix2x2::HermATimesB(temp2, temp, temp);
				//temp2[0] = 1.0; temp2[1] = 0.0; temp[2] = 0.0; temp[3] = 1.0;
				//Matrix2x2::MultiplyAdd(measWeightIter, temp2, weightScalar);
				Matrix2x2::MultiplyAdd(measWeightIter, temp2, weightScalar);
			}
			
			measFluxIter += 4;
			measWeightIter += 4;
			beamWeightPtr += 4;
		}
	}
}
	