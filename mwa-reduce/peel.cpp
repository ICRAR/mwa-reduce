#include "calibrationmethod.h"
#include "model.h"
#include "banddata.h"
#include "solutionfile.h"
#include "beamevaluator.h"
#include "matrix2x2.h"
#include "mspredicter.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <boost/thread/thread.hpp>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <queue>
#include <complex>

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
		if((iters >= data.nIter || !std::isfinite(limit)) && !(*(data.calMethods))[taskIndex]->OnlySolveRotation())
		{
			std::cout << "Recalculating channel " << taskIndex << " (precision=" << limit << ").\n";
			(*(data.calMethods))[taskIndex]->InitSolutionsToUnity();
			iters = data.nIter;
			limit = data.limit;
			(*(data.calMethods))[taskIndex]->Execute(limit, iters);
			if(iters >= data.nIter || !std::isfinite(limit))
			{
				std::cout << "Channel " << taskIndex << " did not converge, setting gains to NaN.\n";
				(*(data.calMethods))[taskIndex]->InitSolutionsToNaN();
			}
			else {
				lastTaskIndex = taskIndex;
			}
		}
		else {
			lastTaskIndex = taskIndex;
		}
		lock.lock();
	}
}

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: peel [-datacolumn <column>] [-beam-on-source] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist>] [-l <precision>] [-i <niter>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-applybeam] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n";
	} else {
		int argi = 1;
		bool
			savePlotFiles = false, saveFaradayPlotFiles = false, saveCrossTermsPlotFile = false, beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, modelFile, rhsSolutionFile;
		std::string dataColumnName = "DATA";
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
			else if(strcmp(argv[argi], "-pf") == 0)
			{
				saveFaradayPlotFiles = true;
				plotFaradayFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-px") == 0)
			{
				saveCrossTermsPlotFile = true;
				crossTermsPlotFile = argv[argi+1];
				argi += 2;
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
			else if(strcmp(argv[argi], "-applybeam") == 0)
			{
				applyBeam = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-beam-on-source") == 0)
			{
				beamOnSource = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-scalar") == 0)
			{
				onlyScalar = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-diag") == 0)
			{
				onlyDiag = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-rhs") == 0)
			{
				rhsSolutionFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-rotation") == 0)
			{
				onlyRotation = true;
				argi++;
			}
			else if(strcmp(argv[argi], "-datacolumn") == 0)
			{
				++argi;
				dataColumnName = argv[argi];
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
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
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROArrayColumn<complex_t> dataColumn(ms, dataColumnName);
		casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
		casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		std::unique_ptr<casa::ROArrayColumn<complex_t>> modelColumn;
		
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		if(polarizationCount != 4)
			throw std::runtime_error("Pol count in MS != 4");
		
		std::cout << "DONE\nCounting timesteps... " << std::flush;
		double time = -1.0;
		std::vector<size_t> timestepRows;
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			if(timeColumn(rowIndex) != time)
			{
				timestepRows.push_back(rowIndex);
				time = timeColumn(rowIndex);
			}
		}
		size_t timestepCount = timestepRows.size();
		timestepRows.push_back(ms.nrow());
		std::cout << "DONE (" << timestepCount << ")\n";
		
		size_t timestepsPerSolution = 1;
		
		size_t passCount = (timestepCount + timestepsPerSolution - 1) / timestepsPerSolution;
		
		std::unique_ptr<Model> model;
		if(!modelFile.empty()) {
				std::cout << "Reading model... " << std::flush;
				model.reset(new Model(modelFile.c_str()));
				std::cout << "DONE\n";
		}

		for(size_t ant=0; ant!=antennaCount; ++ant)
		{
			std::ostringstream antFilename;
			antFilename << "peel-sol-ant" << ant << ".txt";
			std::ofstream(antFilename.str().c_str());
		}
			
		SolutionFile solutionFile;
		solutionFile.SetAntennaCount(antennaCount);
		solutionFile.SetChannelCount(channelCount);
		solutionFile.SetPolarizationCount(4);
		solutionFile.OpenForWriting(outName);

		for(size_t pass=0; pass!=passCount; ++pass)
		{
			size_t
				startTimestep = timestepCount * pass / passCount,
				endTimestep = timestepCount * (pass+1) / passCount,
				timestepsInPass = endTimestep - startTimestep,
				startRow = timestepRows[startTimestep],
				endRow = timestepRows[endTimestep];
			std::vector<CalibrationMethod*> calMethods(channelCount);
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				calMethods[ch] = new CalibrationMethod(1, antennaCount, timestepsInPass);
				calMethods[ch]->SetOnlySolveScalar(onlyScalar);
				calMethods[ch]->SetOnlySolveDiag(onlyDiag);
				calMethods[ch]->SetOnlySolveRotation(onlyRotation);
			}
			std::unique_ptr<MSPredicter> predicter;
			std::unique_ptr<BeamEvaluator> beamEvaluator;
			std::vector<std::complex<double>> beamValues;
			if(modelFile.empty()) {
				std::cout << "Reading data and model column... " << std::flush;
				modelColumn.reset(new casa::ROArrayColumn<complex_t>(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
			}
			else {
				if(beamOnSource || applyBeam)
				{
					beamEvaluator.reset(new BeamEvaluator(ms));
				}
				if(beamOnSource)
				{
					if(model->SourceCount() != 1)
						std::cout << "Warning: To correct for the beam, there should be exactly one source in the model";
					const ModelSource& source = model->Source(0);
					beamValues.resize(channelCount*4);
					double beamSum[4] = {0.0, 0.0, 0.0, 0.0};
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						double frequency = bandData.ChannelFrequency(ch);
						beamEvaluator->EvaluateApparentToAbsGain(source.PosRA(), source.PosDec(), frequency, &beamValues[ch*4]);
						for(size_t p=0; p!=4; ++p)
							beamSum[p] += std::abs(beamValues[ch*4+p]);
					}
				}
				
				predicter.reset(new MSPredicter(ms, *model));
				predicter->SetApplyBeam(applyBeam);
				predicter->SetStartRow(startRow);
				predicter->SetEndRow(endRow);
			}
			
			std::vector<std::complex<double> > modelValues(4 * channelCount);
			casa::Array<complex_t> data(dataShape), modelData(dataShape);
			casa::Array<float> weights(dataShape);
			casa::Array<bool> flags(dataShape);
			size_t timeIndex = 0;
			time = timeColumn(startRow);
			size_t selectedCount = 0, notSelected = 0;
			MSPredicter::RowData rowData;
			
			predicter->Start(false);
			while(predicter->GetNextRow(rowData))
			{
				size_t rowIndex = rowData.rowIndex;
				boost::mutex::scoped_lock lock(predicter->IOMutex());
				if(timeColumn(rowIndex) != time)
				{
					++timeIndex;
					time = timeColumn(rowIndex);
				}
				// Cross correlation?
				size_t antenna1 = rowData.a1, antenna2 = rowData.a2;
				if(antenna1 == antenna2)
				{
					lock.unlock();
				}
				else {
					dataColumn.get(rowIndex, data);
					weightColumn.get(rowIndex, weights);
					flagColumn.get(rowIndex, flags);
					if(modelFile.empty())
						modelColumn->get(rowIndex, modelData);
					lock.unlock();
					
					std::complex<float> *dataPtr = data.cbegin();
					float *weightsPtr = weights.cbegin();
					bool *flagPtr = flags.cbegin();
				
					double u = rowData.u;
					double v = rowData.v;
					double w = rowData.w;
					
					bool selected = true;
					if(u*u + v*v + w*w < minUVW*minUVW)
						selected = false;
					if(selected)
						selectedCount++;
					else
						notSelected++;
				
					if(modelFile.empty())
					{
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
						for(size_t ch = 0; ch!=channelCount; ++ch)
						{
							size_t chIndex = ch * 4;
							for(size_t p=0; p!=4; ++p)
							{
								modelValues[chIndex+p] = rowData.modelData[chIndex+p];
								if(flagPtr[chIndex+p] || !selected) weightsPtr[chIndex+p] = 0.0;
							}
							if(beamOnSource)
							{
								std::complex<double>
									tempResult[4],
									doubleData[4] = {dataPtr[chIndex],dataPtr[chIndex+1],dataPtr[chIndex+2],dataPtr[chIndex+3]};
								Matrix2x2::ATimesB(tempResult, &beamValues[ch*4], doubleData);
								Matrix2x2::ATimesHermB(doubleData, tempResult, &beamValues[ch*4]);
								calMethods[ch]->AddData(doubleData, &weightsPtr[chIndex], &modelValues[chIndex], antenna1, antenna2, timeIndex);
							}
							else {
								calMethods[ch]->AddData(&dataPtr[chIndex], &weightsPtr[chIndex], &modelValues[chIndex], antenna1, antenna2, timeIndex);
							}
						}
					}
				}
				
				predicter->FinishRow(rowData);
			}
		
			std::queue<size_t> tasks;
			for(size_t ch=0; ch!=channelCount; ++ch)
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
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					std::complex<double> val[4];
					for(size_t p=0; p!=4; ++p)
						val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
					Matrix2x2::Invert(val);
					
					for(size_t p=0; p!=4; ++p)
					{
						solutionFile.WriteSolution(val[p], ant, ch, p);
					}
				}
			}
			
			double refPhaseXX = 0.0, refPhaseYY = 0.0;
			for(size_t ant=0; ant!=antennaCount; ++ant)
			{
				std::ostringstream antFilename;
				antFilename << "peel-sol-ant" << ant << ".txt";
				std::ofstream antFile(antFilename.str().c_str(), std::ios_base::app);
				double sumPhaseXX = 0.0, sumPhaseYY = 0.0, sumGainXX = 0.0, sumGainYY = 0.0;
				size_t sumCount = 0;
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					std::complex<double> val[4];
					for(size_t p=0; p!=4; ++p)
						val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
					Matrix2x2::Invert(val);
					if(std::isfinite(val[0].real()) && std::isfinite(val[3].real()) &&
						std::isfinite(val[0].imag()) && std::isfinite(val[3].imag()))
					{
						sumGainXX += std::abs(val[0]);
						sumGainYY += std::abs(val[3]); 
						sumPhaseXX += std::arg(val[0]);
						sumPhaseYY += std::arg(val[3]);
						++sumCount;
					}
				}
				if(ant==0)
				{
					refPhaseXX = sumPhaseXX/sumCount;
					refPhaseYY = sumPhaseYY/sumCount;
				}
				sumPhaseXX = (sumPhaseXX/sumCount - refPhaseXX) * (180.0 / M_PI);
				sumPhaseYY = (sumPhaseYY/sumCount - refPhaseYY) * (180.0 / M_PI);
				
				antFile << startTimestep
					<< '\t' << (sumGainXX/sumCount) << '\t' << sumPhaseXX << '\t'
					<< '\t' << (sumGainYY/sumCount) << '\t' << sumPhaseYY << '\n';
				if(ant == 1)
					std::cout << startTimestep
						<< '\t' << (sumGainXX/sumCount) << '\t' << sumPhaseXX << '\t'
						<< '\t' << (sumGainYY/sumCount) << '\t' << sumPhaseYY << '\n';
			}
				
			if(savePlotFiles)
			{
				std::ofstream phasePlotStream(plotPhaseFile.c_str()), gainPlotStream(plotGainFile.c_str());
				phasePlotStream << antennaCount << ' ' << channelCount << " 4\n";
				gainPlotStream << antennaCount << ' ' << channelCount << " 4\n";
				
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					phasePlotStream << ch << '\t';
					gainPlotStream << ch << '\t';
					
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
			
			if(saveFaradayPlotFiles)
			{
				std::ofstream faradayPlotStream(plotFaradayFile.c_str());
				
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					faradayPlotStream << ch << '\t';
					
					for(size_t ant=0; ant!=antennaCount; ++ant)
					{
						std::complex<double> val[4];
						for(size_t p=0; p!=4; ++p)
							val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
				
						faradayPlotStream << '\t' << -Matrix2x2::RotationAngle(val);
					}
					faradayPlotStream << '\n';
				}
			}
			
			if(saveCrossTermsPlotFile)
			{
				std::ofstream crossTermPlotStream(crossTermsPlotFile.c_str());
				
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					crossTermPlotStream << ch << '\t';
					
					for(size_t ant=0; ant!=antennaCount; ++ant)
					{
						std::complex<double> val[4];
						for(size_t p=0; p!=4; ++p)
							val[p] = calMethods[ch]->JonesSolution(ant, 0, p);
						Matrix2x2::Invert(val);
						double totalPower = std::abs(val[0]) + std::abs(val[1]) + std::abs(val[2]) + std::abs(val[3]);
						crossTermPlotStream << '\t' << (std::abs(val[1]) + std::abs(val[2]))*100.0/totalPower;
					}
					crossTermPlotStream << '\n';
				}
			}
			
			for(size_t ch=0; ch!=channelCount; ++ch)
				delete calMethods[ch];
		}
	}
}
