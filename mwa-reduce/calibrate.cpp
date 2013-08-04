#include "calibrationmethod.h"
#include "model.h"
#include "banddata.h"
#include "solutionfile.h"
#include "beamevaluator.h"
#include "matrix2x2.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <boost/thread/thread.hpp>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <queue>

struct ThreadData
{
	ThreadData() { }
	
	boost::mutex *mutex;
	std::queue<size_t> *tasks;
	std::vector<CalibrationMethod*> *calMethods;
	double limit;
	size_t nIter;
};

void ThreadFunction(ThreadData data)
{
	std::cout << "Thread signing up for duty!\n";
	boost::mutex::scoped_lock lock(*data.mutex);
	size_t lastTaskIndex = data.tasks->front();
	while(!data.tasks->empty()) {
		size_t taskIndex = data.tasks->front();
		data.tasks->pop();
		lock.unlock();
		if(lastTaskIndex != taskIndex)
			(*(data.calMethods))[taskIndex]->InitSolutions(*(*(data.calMethods))[lastTaskIndex]);
		size_t iters = data.nIter;
		double limit = data.limit;
		(*(data.calMethods))[taskIndex]->Execute(limit, iters);
		if(iters >= data.nIter)
		{
			std::cout << "Recalculating channel " << taskIndex << " (precision=" << limit << ").\n";
			(*(data.calMethods))[taskIndex]->InitSolutionsToUnity();
			iters = data.nIter;
			limit = data.limit;
			(*(data.calMethods))[taskIndex]->Execute(data.limit, iters);
		}
		lastTaskIndex = taskIndex;
		lock.lock();
		if(taskIndex<=16)
		{
			std::cout << "Current value of Jones matrix for ant 1, ch " << taskIndex << ":\n"
			<< CalibrationMethod::MatrixToString(& (*(data.calMethods))[taskIndex]->JonesSolution(1, 0, 0));
		}
	
		std::cout << "Finished calibrating channel " << taskIndex << " in " << iters << " iterations, precision=" << limit << ".\n";
	}
	std::cout << "Thread done.\n";
}

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: calibrate [-p <phases.txt> <gains.txt>] [-minuv <min uvw dist>] [-l <precision>] [-i <niter>] [-m <model>] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n"
			<< "Option -a will average over frequency before fitting, nr should specify the amount\n"
			<< "of desired channels.\n";
	} else {
		int argi = 1;
		bool savePlotFiles = false;
		std::string plotPhaseFile, plotGainFile, modelFile;
		size_t niter = 25;
		double limit = 0.0001, minUVW = 0.0;
		
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-p") == 0)
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 3;
			}
			else if(strcmp(argv[argi], "-i") == 0)
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-l") == 0)
			{
				limit = atof(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-m") == 0)
			{
				modelFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-minuv") == 0)
			{
				minUVW = atof(argv[argi+1]);
				argi += 2;
			}
			else throw std::runtime_error("Invalid parameter");
		}
		
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");
		
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		casa::MeasurementSet ms(msName);
		
		std::cout << "Reading meta data... " << std::flush;
		
		/**
		 * Read some meta data from the measurement set
		 */
		casa::MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::MSField fieldTable = ms.field();
		casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		casa::Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
		casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
		casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
		std::unique_ptr<casa::ROArrayColumn<complex_t>> modelColumn;
		
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		if(polarizationCount != 4)
			throw std::runtime_error("Pol count in MS != 4");
		
		std::cout << "DONE\nCounting timesteps... " << std::flush;
		double time = -1.0;
		size_t timestepCount = 0;
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			if(timeColumn(rowIndex) != time)
			{
				++timestepCount;
				time = timeColumn(rowIndex);
			}
		}
		std::cout << "DONE (" << timestepCount << ")\n";

		long int
			pageCount = sysconf(_SC_PHYS_PAGES),
			pageSize = sysconf(_SC_PAGE_SIZE);
		int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
		double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);
		std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";

		size_t nBaselines = antennaCount * (antennaCount-1) / 2;
		size_t samplesPerChannel = nBaselines * timestepCount * 4;
		// 2 for complex data, 2 for complex model, 1 for weights
		double memPerChannel = samplesPerChannel * 5 * sizeof(double);
		std::cout << "One channel takes " << round(memPerChannel*10.0/(1024*1024))/10.0 << " MB of mem.\n";
		size_t channelsPerPass = memSize / memPerChannel;
		if(channelsPerPass > channelCount)
			channelsPerPass = channelCount;
		if(channelsPerPass == 0) {
			std::cout << "WARNING: NOT ENOUGH MEMORY FOR EVEN ONE CHANNEL, expect very bad performance.\n";
			channelsPerPass = 1;
		}
		size_t passCount = channelCount / channelsPerPass;
		std::cout << "Number of channels that fit in memory: " << channelsPerPass << '\n';
		
		std::unique_ptr<Model> model;
		if(!modelFile.empty()) {
				std::cout << "Reading model... " << std::flush;
				model.reset(new Model(modelFile.c_str()));
				std::cout << "DONE\n";
		}

		SolutionFile solutionFile;
		solutionFile.SetAntennaCount(antennaCount);
		solutionFile.SetChannelCount(channelCount);
		solutionFile.SetPolarizationCount(4);
		solutionFile.OpenForWriting(outName);

		for(size_t pass=0; pass!=passCount; ++pass) {
			size_t
				startChannel = (channelCount * pass) / passCount,
				endChannel = (channelCount * (pass+1)) / passCount,
				partChannelCount = endChannel - startChannel;

			BandData partBandData(bandData, startChannel, endChannel);

			std::vector<CalibrationMethod*> calMethods(partChannelCount);
			for(size_t ch=0; ch!=partChannelCount; ++ch)
				calMethods[ch] = new CalibrationMethod(1, antennaCount, timestepCount);
			std::unique_ptr<Predicter> predicter;
			std::unique_ptr<BeamEvaluator> beamEvaluator;
			if(modelFile.empty()) {
				std::cout << "Reading data and model column... " << std::flush;
				modelColumn.reset(new casa::ROArrayColumn<complex_t>(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
			}
			else {
				beamEvaluator.reset(new BeamEvaluator(ms));
				predicter.reset(new Predicter(phaseCentreRA, phaseCentreDec, partBandData.LowestFrequency(), partBandData.HighestFrequency(), partChannelCount, true));
				predicter->Initialize(*model, &*beamEvaluator);
				predicter->ReportSources(*model);
				std::cout << "Reading data & predicting model... " << std::flush;
			}
		
			std::vector<std::complex<double> > modelValues(4 * channelCount);
			casa::Array<complex_t> data(dataShape), modelData(dataShape);
			casa::Array<float> weights(dataShape);
			casa::Array<bool> flags(dataShape);
			size_t timeIndex = 0;
			time = timeColumn(0);
			size_t selectedCount = 0, notSelected = 0;
			for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
			{
				if(timeColumn(rowIndex) != time)
				{
					++timeIndex;
					time = timeColumn(rowIndex);
					std::cout << '.' << std::flush;
				}
				// Cross correlation?
				size_t antenna1 = ant1Column.get(rowIndex), antenna2 = ant2Column.get(rowIndex);
				if(antenna1 != antenna2)
			  {
					dataColumn.get(rowIndex, data);
					weightColumn.get(rowIndex, weights);
					flagColumn.get(rowIndex, flags);
					std::complex<float> *dataPtr = data.cbegin();
					float *weightsPtr = weights.cbegin();
					bool *flagPtr = flags.cbegin();
				
					casa::Array<double> uvwArray = uvwColumn(rowIndex);
					casa::Array<double>::const_iterator i = uvwArray.begin();
					double u = *i; ++i;
					double v = *i; ++i;
					double w = *i;
					
					bool selected = true;
					if(u*u + v*v + w*w < minUVW*minUVW)
						selected = false;
					if(selected)
						selectedCount++;
					else
						notSelected++;
				
					if(modelFile.empty())
				  {
						modelColumn->get(rowIndex, modelData);
						std::complex<float> *modelDataPtr = modelData.cbegin();
						for(size_t ch = 0; ch!=channelCount; ++ch)
					  {
							for(size_t p=0; p!=4; ++p)
						  {
								modelValues[ch*4+p] = *modelDataPtr;
								++modelDataPtr;
							}
						}
					}
					else {
						for(size_t ch = 0; ch!=partChannelCount; ++ch)
					  {
							double lambda = partBandData.ChannelWavelength(ch);
							size_t chIndex = (ch + startChannel) * 4;
							for(size_t p=0; p!=4; ++p)
						  {
								if(flagPtr[chIndex+p] || !selected) weightsPtr[chIndex+p] = 0.0;
								std::complex<double> pVal = predicter->Predict(*model, u/lambda,  v/lambda, w/lambda, ch, p);
								modelValues[chIndex+p] = pVal;
						  }
						  
						  /*modelValues[chIndex+0] = 1.0;
						  modelValues[chIndex+1] = 0.0;
							modelValues[chIndex+2] = 0.0;
							modelValues[chIndex+3] = 1.0;
						  dataPtr[chIndex+0] = 0.1;
							dataPtr[chIndex+1] = 0.0;
							dataPtr[chIndex+2] = 0.0;
							dataPtr[chIndex+3] = 0.1;*/
						  
							calMethods[ch]->AddData(&dataPtr[chIndex], &weightsPtr[chIndex], &modelValues[chIndex], antenna1, antenna2, timeIndex);
						}
					}					
				}
			}
			std::cout << "DONE (" << selectedCount<< "/" << (selectedCount+notSelected) << " rows selected)\nCalibrating...\n";
		
			std::queue<size_t> tasks;
			for(size_t ch=0; ch!=partChannelCount; ++ch)
				tasks.push(ch);
			size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
			boost::thread_group threadGroup;
			boost::mutex mutex;
			for(size_t i=0; i!=cpuCount; ++i)
			{
				ThreadData threadData;
				threadData.mutex = &mutex;
				threadData.tasks = &tasks;
				threadData.calMethods = &calMethods;
				threadData.limit = limit;
				threadData.nIter = niter;
				threadGroup.add_thread(new boost::thread(ThreadFunction, threadData));
			}
			threadGroup.join_all();

			for(size_t ant=0; ant!=antennaCount; ++ant)
		  {
				for(size_t ch=0; ch!=partChannelCount; ++ch)
			  {
					std::complex<double> val[4];
					for(size_t p=0; p!=4; ++p)
						val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
					Matrix2x2::Invert(val);
					
					for(size_t p=0; p!=4; ++p)
				  {
						solutionFile.WriteSolution(val[p], ant, ch+startChannel, p);
					}
				}
			}
			
			if(savePlotFiles)
		  {
				std::ofstream phasePlotStream(plotPhaseFile.c_str()), gainPlotStream(plotGainFile.c_str());
				phasePlotStream << antennaCount << ' ' << partChannelCount << " 4\n";
				gainPlotStream << antennaCount << ' ' << partChannelCount << " 4\n";
				
				for(size_t ch=0; ch!=partChannelCount; ++ch)
			  {
					phasePlotStream << (ch+startChannel) << '\t';
					gainPlotStream << (ch+startChannel) << '\t';
					
					for(size_t p=0; p!=4; ++p)
				  {
						for(size_t ant=0; ant!=antennaCount; ++ant)
					  {
							std::complex<double> val[4];
							for(size_t p2=0; p2!=4; ++p2)
								val[p2] = calMethods[ch]->JonesSolution(ant, 0, p2);
							Matrix2x2::Invert(val);
					
							double s1, s2;
							Matrix2x2::SingularValues(val, s1, s2);
							switch(p)
							{
							case 0: gainPlotStream << '\t' << s1; break;
							case 1: case 2: gainPlotStream << '\t' << 0.0; break;
							case 3: gainPlotStream << '\t' << s2;
							}
							phasePlotStream << '\t' << std::arg(val[p]);
						}
					}
					phasePlotStream << '\n';
					gainPlotStream << '\n';
				}
			}
			
			for(size_t ch=0; ch!=partChannelCount; ++ch)
				delete calMethods[ch];
		}
	}
}
