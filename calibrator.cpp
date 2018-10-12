#include "calibrator.h"

#include "calibrationmethod.h"
#include "banddata.h"
#include "matrix2x2.h"
#include "parallelfor.h"

#include "beamevaluator.h"
#include "mspredicter.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <complex>
#include <fstream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <thread>

Calibrator::Calibrator(casacore::MeasurementSet& ms, size_t threadCount) :
	_ms(ms),
	_dataColumnName("DATA"),
	_minAccuracy(CalibrationMethod::DefaultMinAccuracy()),
	_stoppingAccuracy(CalibrationMethod::DefaultStoppingAccuracy()),
	_nIter(1000),
	_solutionInterval(0),
	_threadCount(threadCount),
	_onlyScalar(false),
	_onlyDiag(false),
	_onlyRotation(false),
	_applyBeam(false),
	_minUVW(0.0),
	_maxUVW(100000000.0),
	_savePlotFiles(false),
	_saveFaradayPlotFiles(false),
	_saveCrossTermsPlotFile(false),
	_verbose(false),
	_selection(MSSelection::Everything()),
	_followAntenna(0)
{
}

void Calibrator::Perform()
{
	if(_verbose)
		std::cout << "Reading meta data... " << std::flush;
	
	/**
		* Read some meta data from the measurement set
		*/
	casacore::MSAntenna aTable = _ms.antenna();
	size_t antennaCount = aTable.nrow();
	
	BandData bandData(_ms.spectralWindow());
	size_t channelCount = bandData.ChannelCount();
	if(channelCount == 0) throw std::runtime_error("No channels in set");
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	typedef float num_t;
	typedef std::complex<num_t> complex_t;
	casacore::ROScalarColumn<double> timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
	casacore::ROArrayColumn<complex_t> dataColumn(_ms, _dataColumnName);
	casacore::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM));
	casacore::ROArrayColumn<bool> flagColumn(_ms, _ms.columnName(casacore::MSMainEnums::FLAG));
	
	casacore::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	
	if(polarizationCount != 4)
		throw std::runtime_error("Pol count in MS != 4");
	
	if(channelCount % _solutionChannels != 0)
	{
		std::cout
			<< "WARNING: You've requested to solve for " << _solutionChannels << " channels\n"
			<< "per solution. However, the total number of channels (" << channelCount << ") is\n"
			<< "not divisable by this. Therefore, some solutions will have different number of channels.\n";
	}
	size_t chBlockCount = channelCount / _solutionChannels;
	if(chBlockCount == 0)
		throw std::runtime_error("Invalid value for solution-channels");
	
	if(_verbose)
	 std::cout << "DONE\nCounting timesteps... " << std::flush;
	double time = -1.0;
	std::vector<size_t> timestepRows, selectedRows;
	for(size_t rowIndex=0;rowIndex!=_ms.nrow();++rowIndex)
	{
		if(timeColumn(rowIndex) != time)
		{
			if(_selection.IsTimeSelected(timestepRows.size()))
			{
				selectedRows.push_back(rowIndex);
			}
			else if(_selection.HasInterval() && timestepRows.size() == _selection.IntervalEnd())
			{
				selectedRows.push_back(rowIndex);
			}
			timestepRows.push_back(rowIndex);
			time = timeColumn(rowIndex);
		}
	}
	if(_selection.IsTimeSelected(timestepRows.size()-1))
		selectedRows.push_back(_ms.nrow());
	size_t selTimestepCount = selectedRows.size()-1;
	size_t intervalCount = (_solutionInterval!=0) ? (selTimestepCount + _solutionInterval - 1) / _solutionInterval : 1;
	if(_verbose)
	 std::cout << "DONE (" << timestepRows.size() << " timesteps, " << selTimestepCount << " selected, " << intervalCount << " intervals)\n";

	if(!_modelFilename.empty()) {
		if(_verbose)
			std::cout << "Reading model... " << std::flush;
		_model = Model(_modelFilename.c_str());
		if(_verbose)
			std::cout << "DONE\n";
	}

	_solutionFile.SetAntennaCount(antennaCount);
	_solutionFile.SetChannelCount(chBlockCount);
	_solutionFile.SetIntervalCount(intervalCount);
	_solutionFile.SetPolarizationCount(4);
	if(_solutionFilename.empty())
		_solutionFile.OpenInMemory();
	else
		_solutionFile.OpenForWriting(_solutionFilename.c_str());

	long int
		pageCount = sysconf(_SC_PHYS_PAGES),
		pageSize = sysconf(_SC_PAGE_SIZE);
	int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);
	// obey absmem limits
	if (_absmem > 0.0)
	{
		memSizeInGB = _absmem;
		memSize = _absmem * 1024.0 * 1024.0 * 1024.0;
	}
	size_t nBaselines = antennaCount * (antennaCount-1) / 2;
	
	for(size_t intervalIndex=0; intervalIndex!=intervalCount; ++intervalIndex)
	{
		std::cout << " >>> INTERVAL " << (intervalIndex+1) << " / " << intervalCount << " <<<\n";
		size_t
			intervalTimestepStart = (intervalIndex*selTimestepCount) / intervalCount,
			intervalTimestepEnd = ((intervalIndex+1)*selTimestepCount) / intervalCount,
			intervalRowStart = selectedRows[intervalTimestepStart],
			intervalRowEnd = selectedRows[intervalTimestepEnd],
			timestepsInInterval = intervalTimestepEnd - intervalTimestepStart;
		size_t samplesPerChannel = nBaselines * timestepsInInterval * 4;
		// 2 for complex data, 2 for complex model, 1 for weights
		double memPerChannelBlock = samplesPerChannel * 5 * sizeof(double) * _solutionChannels;
		if(_verbose)
		{
			std::cout << "Will use " << _threadCount << " cores.\n";
			std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";
			std::cout << "One channel block of " << _solutionChannels << " channels takes " << round(memPerChannelBlock*10.0/(1024*1024))/10.0 << " MB of mem.\n";
		}
		size_t chBlocksPerPass = memSize / memPerChannelBlock;
		if(chBlocksPerPass > chBlockCount)
			chBlocksPerPass = chBlockCount;
		if(chBlocksPerPass == 0)
		{
			if(_verbose)
				std::cout << "WARNING: NOT ENOUGH MEMORY FOR EVEN ONE CHANNEL BLOCK, expect very bad performance.\n";
			chBlocksPerPass = 1;
		}
		size_t passCount = (chBlockCount + chBlocksPerPass - 1) / chBlocksPerPass;
		if(_verbose)
			std::cout << "Number of channel that fit in memory: " << chBlocksPerPass << " blocks (" << chBlocksPerPass*_solutionChannels << " channels, " << passCount << " passes)\n";
		
		for(size_t pass=0; pass!=passCount; ++pass) {
			size_t
				startChBlock = (chBlockCount * pass) / passCount,
				endChBlock = (chBlockCount * (pass+1)) / passCount,
				partChBlockCount = endChBlock - startChBlock;

			_calMethods.resize(partChBlockCount);
			for(size_t cb = 0; cb!=partChBlockCount; ++cb)
			{
				const size_t
					cbChStart = (cb+startChBlock)*channelCount/chBlockCount,
					cbChEnd = (cb+1+startChBlock)*channelCount/chBlockCount,
					cbNCh = cbChEnd - cbChStart;
				std::unique_ptr<CalibrationMethod>& method = _calMethods[cb];
				method.reset(new CalibrationMethod(cbNCh, antennaCount, timestepsInInterval));
				method->SetOnlySolveScalar(_onlyScalar);
				method->SetOnlySolveDiag(_onlyDiag);
				method->SetOnlySolveRotation(_onlyRotation);
			}
			std::unique_ptr<MSPredicter> predicter;
			std::unique_ptr<BeamEvaluator> beamEvaluator;
			std::vector<std::complex<double>> beamValues;
			if(_model.Empty()) {
				std::cout << "Reading data and model column...\n";
				predicter.reset(new MSPredicter(_ms, _threadCount));
				predicter->SetStartRow(intervalRowStart);
				predicter->SetEndRow(intervalRowEnd);
			}
			else {
				if(_applyBeam)
				{
					beamEvaluator.reset(new BeamEvaluator(_ms, _verbose, _mwaPath));
				}
				predicter.reset(new MSPredicter(_ms, _threadCount, _model));
				predicter->SetStartRow(intervalRowStart);
				predicter->SetEndRow(intervalRowEnd);
				predicter->SetApplyBeam(_applyBeam);
				predicter->SetMWAPath(_mwaPath);
				if(_verbose)
					std::cout << "Reading data & predicting model...\n";
			}
			
			std::vector<std::complex<double> > modelValues(4 * channelCount);
			casacore::Array<complex_t> data(dataShape);
			casacore::Array<float> weights(dataShape);
			casacore::Array<bool> flags(dataShape);
			time = timeColumn(intervalRowStart);
			size_t selectedCount = 0, notSelected = 0, previousTime = 0;
			MSPredicter::RowData rowData;
			
			predicter->Start(_verbose);
			while(predicter->GetNextRow(rowData))
			{
				size_t rowIndex = rowData.rowIndex;
				// Cross correlation?
				size_t antenna1 = rowData.a1, antenna2 = rowData.a2;
				if(previousTime < rowData.timeIndex)
				{
					previousTime = rowData.timeIndex;
					if(_verbose)
						std::cout << '.' << std::flush;
				}
				if(antenna1 != antenna2)
				{
					std::unique_lock<std::mutex> lock(predicter->IOMutex());
					dataColumn.get(rowIndex, data);
					weightColumn.get(rowIndex, weights);
					flagColumn.get(rowIndex, flags);
					lock.unlock();
					
					std::complex<float> *dataPtr = data.cbegin();
					float *weightsPtr = weights.cbegin();
					bool *flagPtr = flags.cbegin();
				
					double u = rowData.u;
					double v = rowData.v;
					double w = rowData.w;
					
					bool selected = true;
					if((u*u + v*v + w*w < _minUVW*_minUVW) || (u*u + v*v + w*w > _maxUVW*_maxUVW))
						selected = false;
					if(selected)
						selectedCount++;
					else
						notSelected++;
				
					for(size_t cb = 0; cb!=partChBlockCount; ++cb)
					{
						size_t cbStartCh = (cb + startChBlock)*channelCount/chBlockCount;
						size_t cbEndCh = (cb + 1 + startChBlock)*channelCount/chBlockCount;
						size_t cbChCount = cbEndCh - cbStartCh;
						for(size_t ch = 0; ch!=cbChCount; ++ch)
						{
							size_t chIndex = (ch + cbStartCh) * 4;
							for(size_t p=0; p!=4; ++p)
							{
								modelValues[chIndex+p] = rowData.modelData[chIndex+p];
								if(flagPtr[chIndex+p] || !selected) weightsPtr[chIndex+p] = 0.0;
							}
							_calMethods[cb]->AddData(&dataPtr[chIndex], &weightsPtr[chIndex], &modelValues[chIndex], antenna1, antenna2, rowData.timeIndex);
						}
					}
				}
				
				predicter->FinishRow(rowData);
			}
			if(_verbose)
				std::cout << "DONE (" << selectedCount<< "/" << (selectedCount+notSelected) << " rows selected)\nCalibrating...\n";
		
			ParallelFor<size_t> loop(_threadCount);
			_threadData.assign(loop.NThreads(), ThreadData(chBlockCount));
			loop.Run(0, partChBlockCount, std::bind(&Calibrator::calibrateChannelBlock, this, std::placeholders::_1, std::placeholders::_2));

			// Save solutions
			for(size_t ant=0; ant!=antennaCount; ++ant)
			{
				for(size_t cb=0; cb!=partChBlockCount; ++cb)
				{
					std::complex<double> val[4];
					for(size_t p=0; p!=4; ++p)
						val[p] = _calMethods[cb]->JonesSolution(ant, 0, p);
					Matrix2x2::Invert(val);
					
					for(size_t p=0; p!=4; ++p)
					{
						_solutionFile.WriteSolution(val[p], intervalIndex, ant, cb+startChBlock, p);
					}
				}
			}
			
			if(_savePlotFiles)
			{
				std::ofstream phasePlotStream(_phasePlotFilename.c_str()), gainPlotStream(_gainPlotFilename.c_str());
				phasePlotStream << antennaCount << ' ' << partChBlockCount << " 4\n";
				gainPlotStream << antennaCount << ' ' << partChBlockCount << " 4\n";
				
				for(size_t cb=0; cb!=partChBlockCount; ++cb)
				{
					phasePlotStream << (cb+startChBlock) << '\t';
					gainPlotStream << (cb+startChBlock) << '\t';
					
					for(size_t p=0; p!=4; ++p)
					{
						for(size_t ant=0; ant!=antennaCount; ++ant)
						{
							std::complex<double> val[4];
							for(size_t p2=0; p2!=4; ++p2)
								val[p2] = _calMethods[cb]->JonesSolution(ant, 0, p2);
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
				
				for(size_t cb=0; cb!=partChBlockCount; ++cb)
				{
					faradayPlotStream << (cb+startChBlock) << '\t';
					
					for(size_t ant=0; ant!=antennaCount; ++ant)
					{
						std::complex<double> val[4];
						for(size_t p=0; p!=4; ++p)
							val[p] = _calMethods[cb]->JonesSolution(ant, 0, p);
				
						faradayPlotStream << '\t' << -Matrix2x2::RotationAngle(val);
					}
					faradayPlotStream << '\n';
				}
			}
			
			if(_saveCrossTermsPlotFile)
			{
				std::ofstream crossTermPlotStream(_crossTermsPlotFilename.c_str());
				
				for(size_t cb=0; cb!=partChBlockCount; ++cb)
				{
					crossTermPlotStream << (cb+startChBlock) << '\t';
					
					for(size_t ant=0; ant!=antennaCount; ++ant)
					{
						std::complex<double> val[4];
						for(size_t p=0; p!=4; ++p)
							val[p] = _calMethods[cb]->JonesSolution(ant, 0, p);
						Matrix2x2::Invert(val);
						double totalPower = std::abs(val[0]) + std::abs(val[1]) + std::abs(val[2]) + std::abs(val[3]);
						crossTermPlotStream << '\t' << (std::abs(val[1]) + std::abs(val[2]))*100.0/totalPower;
					}
					crossTermPlotStream << '\n';
				}
			}
			
			_calMethods.clear(); // destructs calibration methods
		}
	}
}

void Calibrator::calibrateChannelBlock(size_t channelBlockIndex, size_t threadIndex)
{
	size_t lastSuccessfulChBlock = _threadData[threadIndex].lastSuccesfulChBlock;
	
	CalibrationMethod& method = *_calMethods[channelBlockIndex];
	
	if(lastSuccessfulChBlock < channelBlockIndex)
		method.InitSolutions(*_calMethods[lastSuccessfulChBlock]);
	size_t iters = _nIter;
	double limit = _stoppingAccuracy;
	method.Execute(limit, iters);
	if((iters >= _nIter || !std::isfinite(limit)) && !method.OnlySolveRotation())
	{
		std::cout << "Recalculating channel " << channelBlockIndex << " (accuracy=" << limit << ").\n";
		method.InitSolutionsToUnity();
		iters = _nIter;
		limit = _stoppingAccuracy;
		method.Execute(limit, iters);

		if((iters >= _nIter && limit > _minAccuracy) || !std::isfinite(limit))
		{
			std::cout << "Channel " << channelBlockIndex << " did not converge (accuracy=" << limit << "), setting gains to NaN.\n";
			method.InitSolutionsToNaN();
		}
		else {
			if(iters >= _nIter && limit > _stoppingAccuracy)
			{
				std::cout << "Channel " << channelBlockIndex << " converged (accuracy=" << limit << ") but did not reach stopping accuracy.\n";
			}
			_threadData[threadIndex].lastSuccesfulChBlock = channelBlockIndex;
		}
	}
	else {
		_threadData[threadIndex].lastSuccesfulChBlock = channelBlockIndex;
	}
	if(channelBlockIndex<=16)
	{
		if(_verbose)
			std::cout << "Current value of Jones matrix for ant " << _followAntenna << ", ch " << channelBlockIndex << ":\n"
		<< CalibrationMethod::MatrixToString(& method.JonesSolution(_followAntenna, 0, 0));
	}

	if(_verbose)
		std::cout << "Finished calibrating channel " << channelBlockIndex << " in " << iters << " iterations, precision=" << limit << ".\n";
}

