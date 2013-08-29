#include "calibrator.h"

#include "calibrationmethod.h"
#include "banddata.h"
#include "beamevaluator.h"
#include "matrix2x2.h"
#include "mspredicter.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <boost/thread/thread.hpp>

#include <fstream>
#include <stdexcept>
#include <memory>
#include <queue>
#include <complex>

Calibrator::Calibrator(casa::MeasurementSet& ms) :
	_ms(ms),
	_dataColumnName("DATA"),
	_limit(0.0001),
	_nIter(100),
	_solutionInterval(0),
	_onlyScalar(false),
	_onlyDiag(false),
	_onlyRotation(false),
	_beamOnSource(false),
	_applyBeam(false),
	_minUVW(0.0),
	_savePlotFiles(false),
	_saveFaradayPlotFiles(false),
	_saveCrossTermsPlotFile(false),
	_verbose(false)
{
}

void Calibrator::Perform()
{
	if(_verbose)
		std::cout << "Reading meta data... " << std::flush;
	
	/**
		* Read some meta data from the measurement set
		*/
	casa::MSAntenna aTable = _ms.antenna();
	size_t antennaCount = aTable.nrow();
	
	BandData bandData(_ms.spectralWindow());
	size_t channelCount = bandData.ChannelCount();
	if(channelCount == 0) throw std::runtime_error("No channels in set");
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	typedef float num_t;
	typedef std::complex<num_t> complex_t;
	casa::ROScalarColumn<double> timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
	casa::ROArrayColumn<complex_t> dataColumn(_ms, _dataColumnName);
	casa::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(_ms, _ms.columnName(casa::MSMainEnums::FLAG));
	std::unique_ptr<casa::ROArrayColumn<complex_t>> modelColumn;
	
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	
	if(polarizationCount != 4)
		throw std::runtime_error("Pol count in MS != 4");
	
	if(_verbose)
	 std::cout << "DONE\nCounting timesteps... " << std::flush;
	double time = -1.0;
	std::vector<size_t> timestepRows;
	for(size_t rowIndex=0;rowIndex!=_ms.nrow();++rowIndex)
	{
		if(timeColumn(rowIndex) != time)
		{
			timestepRows.push_back(rowIndex);
			time = timeColumn(rowIndex);
		}
	}
	size_t timestepCount = timestepRows.size();
	timestepRows.push_back(_ms.nrow());
	size_t intervalCount = (_solutionInterval!=0) ? (timestepCount + _solutionInterval - 1) / _solutionInterval : 1;
	if(_verbose)
	 std::cout << "DONE (" << timestepCount << " timesteps, " << intervalCount << " intervals)\n";

	if(!_modelFilename.empty()) {
		if(_verbose)
			std::cout << "Reading model... " << std::flush;
		_model = Model(_modelFilename.c_str());
		if(_verbose)
			std::cout << "DONE\n";
	}

	_solutionFile.SetAntennaCount(antennaCount);
	_solutionFile.SetChannelCount(channelCount);
	_solutionFile.SetIntervalCount(intervalCount);
	_solutionFile.SetPolarizationCount(4);
	if(_solutionFilename.empty())
		_solutionFile.OpenInMemory();
	else
		_solutionFile.OpenForWriting(_solutionFilename.c_str());

	for(size_t intervalIndex=0; intervalIndex!=intervalCount; ++intervalIndex)
	{
		size_t
			intervalTimestepStart = (intervalIndex*timestepCount) / intervalCount,
			intervalTimestepEnd = ((intervalIndex+1)*timestepCount) / intervalCount,
			intervalRowStart = timestepRows[intervalTimestepStart],
			intervalRowEnd = timestepRows[intervalTimestepEnd],
			timestepsInInterval = intervalTimestepEnd - intervalTimestepStart;
		long int
			pageCount = sysconf(_SC_PHYS_PAGES),
			pageSize = sysconf(_SC_PAGE_SIZE);
		int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
		double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);

		size_t nBaselines = antennaCount * (antennaCount-1) / 2;
		size_t samplesPerChannel = nBaselines * timestepsInInterval * 4;
		// 2 for complex data, 2 for complex model, 1 for weights
		double memPerChannel = samplesPerChannel * 5 * sizeof(double);
		if(_verbose)
		{
			std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";
			std::cout << "One channel takes " << round(memPerChannel*10.0/(1024*1024))/10.0 << " MB of mem.\n";
		}
		size_t channelsPerPass = memSize / memPerChannel;
		if(channelsPerPass > channelCount)
			channelsPerPass = channelCount;
		if(channelsPerPass == 0) {
			if(_verbose)
				std::cout << "WARNING: NOT ENOUGH MEMORY FOR EVEN ONE CHANNEL, expect very bad performance.\n";
			channelsPerPass = 1;
		}
		size_t passCount = (channelCount + channelsPerPass - 1) / channelsPerPass;
		if(_verbose)
			std::cout << "Number of channels that fit in memory: " << channelsPerPass << " (" << passCount << " passes)\n";
		
		for(size_t pass=0; pass!=passCount; ++pass) {
			size_t
				startChannel = (channelCount * pass) / passCount,
				endChannel = (channelCount * (pass+1)) / passCount,
				partChannelCount = endChannel - startChannel;

			BandData partBandData(bandData, startChannel, endChannel);

			std::vector<CalibrationMethod*> calMethods(partChannelCount);
			for(size_t ch=0; ch!=partChannelCount; ++ch)
			{
				calMethods[ch] = new CalibrationMethod(1, antennaCount, timestepsInInterval);
				calMethods[ch]->SetOnlySolveScalar(_onlyScalar);
				calMethods[ch]->SetOnlySolveDiag(_onlyDiag);
				calMethods[ch]->SetOnlySolveRotation(_onlyRotation);
			}
			std::unique_ptr<MSPredicter> predicter;
			std::unique_ptr<BeamEvaluator> beamEvaluator;
			std::vector<std::complex<double>> beamValues;
			if(_model.Empty()) {
				std::cout << "Reading data and model column... " << std::flush;
				modelColumn.reset(new casa::ROArrayColumn<complex_t>(_ms, _ms.columnName(casa::MSMainEnums::MODEL_DATA)));
			}
			else {
				if(_beamOnSource || _applyBeam)
				{
					beamEvaluator.reset(new BeamEvaluator(_ms, _verbose));
				}
				if(_beamOnSource)
				{
					if(_model.SourceCount() != 1)
						std::cout << "Warning: To correct for the beam, there should be exactly one source in the model";
					const ModelSource& source = _model.Source(0);
					std::cout << "Predicting beam... " << std::flush;
					beamValues.resize(partChannelCount*4);
					double beamSum[4] = {0.0, 0.0, 0.0, 0.0};
					for(size_t ch=0; ch!=partChannelCount; ++ch)
					{
						double frequency = partBandData.ChannelFrequency(ch);
						beamEvaluator->EvaluateApparentToAbsGain(source.PosRA(), source.PosDec(), frequency, &beamValues[ch*4]);
						for(size_t p=0; p!=4; ++p)
							beamSum[p] += std::abs(beamValues[ch*4+p]);
					}
					std::cout << "DONE (avg inv beam:";
					for(size_t p=0; p!=4; ++p)
						std::cout << ' ' << beamSum[p]/partChannelCount;
					std::cout << '\n';
				}
				
				predicter.reset(new MSPredicter(_ms, _model));
				predicter->SetStartRow(intervalRowStart);
				predicter->SetEndRow(intervalRowEnd);
				predicter->SetApplyBeam(_applyBeam);
				if(_verbose)
					std::cout << "Reading data & predicting model...\n";
			}
			
			std::vector<std::complex<double> > modelValues(4 * channelCount);
			casa::Array<complex_t> data(dataShape), modelData(dataShape);
			casa::Array<float> weights(dataShape);
			casa::Array<bool> flags(dataShape);
			size_t timeIndex = 0;
			time = timeColumn(intervalRowStart);
			size_t selectedCount = 0, notSelected = 0;
			MSPredicter::RowData rowData;
			
			predicter->Start(_verbose);
			while(predicter->GetNextRow(rowData))
			{
				size_t rowIndex = rowData.rowIndex;
				boost::mutex::scoped_lock lock(predicter->IOMutex());
				if(timeColumn(rowIndex) != time)
				{
					++timeIndex;
					time = timeColumn(rowIndex);
					if(_verbose)
						std::cout << '.' << std::flush;
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
					if(_model.Empty())
						modelColumn->get(rowIndex, modelData);
					lock.unlock();
					
					std::complex<float> *dataPtr = data.cbegin();
					float *weightsPtr = weights.cbegin();
					bool *flagPtr = flags.cbegin();
				
					double u = rowData.u;
					double v = rowData.v;
					double w = rowData.w;
					
					bool selected = true;
					if(u*u + v*v + w*w < _minUVW*_minUVW)
						selected = false;
					if(selected)
						selectedCount++;
					else
						notSelected++;
				
					if(_model.Empty())
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
						for(size_t ch = 0; ch!=partChannelCount; ++ch)
						{
							size_t chIndex = (ch + startChannel) * 4;
							for(size_t p=0; p!=4; ++p)
							{
								modelValues[chIndex+p] = rowData.modelData[chIndex+p];
								if(flagPtr[chIndex+p] || !selected) weightsPtr[chIndex+p] = 0.0;
							}
							if(_beamOnSource)
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
			if(_verbose)
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
				threadData.limit = _limit;
				threadData.nIter = _nIter;
				threadGroup.add_thread(new boost::thread(&Calibrator::threadFunction, this, threadData));
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
						_solutionFile.WriteSolution(val[p], intervalIndex, ant, ch+startChannel, p);
					}
				}
			}
			
			if(_savePlotFiles)
			{
				std::ofstream phasePlotStream(_phasePlotFilename.c_str()), gainPlotStream(_gainPlotFilename.c_str());
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
			
			if(_saveFaradayPlotFiles)
			{
				std::ofstream faradayPlotStream(_faradayPlotFilename.c_str());
				
				for(size_t ch=0; ch!=partChannelCount; ++ch)
				{
					faradayPlotStream << (ch+startChannel) << '\t';
					
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
			
			if(_saveCrossTermsPlotFile)
			{
				std::ofstream crossTermPlotStream(_crossTermsPlotFilename.c_str());
				
				for(size_t ch=0; ch!=partChannelCount; ++ch)
				{
					crossTermPlotStream << (ch+startChannel) << '\t';
					
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
			
			for(size_t ch=0; ch!=partChannelCount; ++ch)
				delete calMethods[ch];
		}
	}

}

void Calibrator::threadFunction(ThreadData data)
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
		if(taskIndex<=16)
		{
			if(_verbose)
				std::cout << "Current value of Jones matrix for ant 1, ch " << taskIndex << ":\n"
			<< CalibrationMethod::MatrixToString(& (*(data.calMethods))[taskIndex]->JonesSolution(1, 0, 0));
		}
	
		if(_verbose)
			std::cout << "Finished calibrating channel " << taskIndex << " in " << iters << " iterations, precision=" << limit << ".\n";
	}
}

