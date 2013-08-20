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

#include <limits>

void SpectrumMaker::measure(const string& filename, const string& solutionsFile)
{
	casa::MeasurementSet ms(filename);
	
	/**
		* Read some meta data from the measurement set
		*/
	BandData bandData(ms.spectralWindow());
	size_t channelCount = bandData.ChannelCount();
	
	casa::MSField fieldTable = ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	
	std::vector<long double>
		measFlux(channelCount * _sources.size() * 4),
		measWeights(channelCount * _sources.size() * 4);
	
	MSPredicter modelPredicter(ms, _subtractedModel, solutionsFile);
	
	std::vector<std::unique_ptr<Predicter>> predicters;
	for(std::vector<ModelSource>::iterator sourceIter=_sources.begin();
			sourceIter!=_sources.end(); ++sourceIter)
	{
		predicters.push_back(std::unique_ptr<Predicter>(
			new Predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount)));
		(*predicters.rbegin())->Initialize(*sourceIter);
	}
	
	const size_t BUFFER_COUNT = 16;
	
	size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	boost::thread_group threadGroup;
	std::vector<lane<ThreadTaskInfo>*> taskLanes(cpuCount);
	for(size_t i=0; i!=cpuCount; ++i)
	{
		// The task lane must issue a "wait" when its size is
		// BUFFER_COUNT-2, because the pushing thread can start
		// overwritting the n-1 th buffer while the popping thread
		// can just get the n-1 th buffer out.
		taskLanes[i] = new lane<ThreadTaskInfo>(BUFFER_COUNT-2);
		threadGroup.add_thread(new boost::thread(&SpectrumMaker::measureThreadFunc, this, &bandData, &*taskLanes[i], channelCount, polarizationCount));
	}
	
	std::vector<std::complex<double>*> dataBuffers(BUFFER_COUNT);
	std::vector<casa::Array<bool>*> flagBuffers(BUFFER_COUNT);
	for(size_t i=0; i!=BUFFER_COUNT; ++i)
	{
		dataBuffers[i] = new std::complex<double>[channelCount * polarizationCount];
		flagBuffers[i] = new casa::Array<bool>(dataShape);
	}
	casa::Array<casa::Complex> dataArray(dataShape);
	
	modelPredicter.Start();
	
	/**
		* Calculate spectra
		*/
	size_t bufferIndex = 0;
	MSPredicter::RowData rowData;
	ProgressBar progress(std::string("Measure spectra in ") + filename);
	while(modelPredicter.GetNextRow(rowData))
	{
		size_t rowIndex = rowData.rowIndex;
		
		boost::mutex::scoped_lock lock(modelPredicter.IOMutex());
		progress.SetProgress(rowIndex, ms.nrow());
		
		// Cross correlation?
		size_t
			a1 = ant1Column.get(rowIndex),
			a2 = ant2Column.get(rowIndex);
		if(a1 != a2)
		{
			std::complex<double> *data = dataBuffers[bufferIndex];
			casa::Array<bool> &flags = *flagBuffers[bufferIndex];
			
			dataColumn.get(rowIndex, dataArray);
			flagColumn.get(rowIndex, flags);
			lock.unlock();
			
			casa::Array<casa::Complex>::const_contiter dataArrayIter = dataArray.cbegin();
			for(size_t i=0; i!=polarizationCount*channelCount; ++i)
			{
				data[i] = std::complex<double>(dataArrayIter->real(), dataArrayIter->imag()) - rowData.modelData[i];
				++dataArrayIter;
			}
			
			for(size_t s=0; s!=_sources.size(); ++s)
			{
				size_t thread = s % cpuCount;
				ThreadTaskInfo task;
				task.flux = &measFlux[s * channelCount * polarizationCount];
				task.fluxWeights = &measWeights[s * channelCount * polarizationCount];
				task.data = data;
				task.flags = flags.cbegin();
				task.u = rowData.u;
				task.v = rowData.v;
				task.w = rowData.w;
				task.sourceIndex = s;
				task.predicter = &*predicters[s];
				task.a1 = a1;
				task.a2 = a2;
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
		delete flagBuffers[i];
	}
	
	progress.SetProgress(ms.nrow(), ms.nrow());
	
	// Apply beam
	BeamEvaluator beamEval(ms);
	for(size_t s=0; s!=_sources.size(); ++s)
	{
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			const ModelSource& src = _sources[s];
			beamEval.ApparentToAbs(src.PosRA(), src.PosDec(), bandData.ChannelFrequency(ch), &measFlux[(s * channelCount + ch) * 4]);
		}
	}
	
	// Add to the total values
	for(size_t s=0; s!=_sources.size(); ++s)
	{
		Spectrum &spectrum = _spectrumPerSource[s];
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double freq = bandData.ChannelFrequency(ch);
			spectrum.AddMeasurement(freq, &measFlux[(s * channelCount + ch) * 4], &measWeights[(s * channelCount + ch) * 4]);
		}
	}
}

void SpectrumMaker::measureThreadFunc(const BandData* bandData, lane< SpectrumMaker::ThreadTaskInfo >* taskLane, size_t channelCount, size_t polarizationCount)
{
	ThreadTaskInfo taskInfo;
	while(taskLane->read(taskInfo))
	{
		std::complex<double>* dataPtr = taskInfo.data;
		bool* flagPtr = taskInfo.flags;
		long double* measFluxIter = taskInfo.flux;
		long double* measWeightIter = taskInfo.fluxWeights;
		Predicter& predicter = *taskInfo.predicter;
		double
			u = taskInfo.u,
			v = taskInfo.v,
			w = taskInfo.w;
		size_t s = taskInfo.sourceIndex;
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double lambda = bandData->ChannelWavelength(ch);
			Predicter::CNumType predicted[4];
			predicter.Predict4(predicted, _sources[s], u/lambda, v/lambda, w/lambda, ch, taskInfo.a1, taskInfo.a2);
			std::complex<double> weight[4]; // TODO determine these
			double visSample[4];
			bool sampleGood = true;
			for(size_t p=0; p!=4; ++p)
			{
				double 
					real = dataPtr->real(),
					imag = dataPtr->imag();
				if(!(*flagPtr) && std::isfinite(real) && std::isfinite(imag))
				{
					// add real(data * conj(predicted))
					visSample[p] = real * predicted[0].real() + imag * predicted[0].imag();
				}
				else {
					sampleGood = false;
				}
				++dataPtr;
				++flagPtr;
			}
			
			if(sampleGood)
			{
				std::complex<double> temp[4];
				Matrix2x2::HermATimesB(temp, weight, visSample);
				Matrix2x2::PlusATimesB(measFluxIter, temp, weight);
				
				Matrix2x2::PlusATimesHermB(measWeightIter, weight, weight);
			}
			
			measFluxIter += 4;
			measWeightIter += 4;
		}
	}
}
	