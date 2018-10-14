#include "peeler.h"

#include "calibrationmethod.h"
#include "solutionfile.h"
#include "beamevaluator.h"
#include "matrix2x2.h"
#include "mspredicter.h"

#include "model/model.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <fstream>
#include <stdexcept>
#include <memory>
#include <complex>

Peeler::Peeler(casacore::MeasurementSet& ms, size_t threadCount) :
	_ms(ms),
	_beamOnSource(false),
	_applyBeam(false),
	_onlyScalar(false),
	_onlyDiag(false),
	_onlyRotation(false),
	_saveSolutionsFiles(false),
	_dataColumnName("DATA"),
	_nIter(CalibrationMethod::DefaultNIter()),
	_minAccuracy(CalibrationMethod::DefaultMinAccuracy()),
	_stoppingAccuracy(CalibrationMethod::DefaultStoppingAccuracy()),
	_minUVW(0.0),
	_solutionInterval(1),
	_threadCount(threadCount)
{
	_ms.reopenRW();
}

void Peeler::Perform()
{
	std::cout << "Reading meta data... " << std::flush;
	
	/**
		* Read some meta data from the measurement set
		*/
	casacore::MSAntenna aTable = _ms.antenna();
	size_t antennaCount = aTable.nrow();
	
	_bandData = BandData(_ms.spectralWindow());
	size_t channelCount = _bandData.ChannelCount();
	if(channelCount == 0) throw std::runtime_error("No channels in set");
	if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	typedef float num_t;
	typedef std::complex<num_t> complex_t;
	casacore::ROScalarColumn<double> timeColumn(_ms, _ms.columnName(casacore::MSMainEnums::TIME));
	casacore::ArrayColumn<complex_t> dataColumn(_ms, _dataColumnName);
	casacore::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM));
	casacore::ArrayColumn<bool> flagColumn(_ms, _ms.columnName(casacore::MSMainEnums::FLAG));
	std::unique_ptr<casacore::ROArrayColumn<complex_t>> modelColumn;
	
	casacore::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	
	if(polarizationCount != 4)
		throw std::runtime_error("Pol count in MS != 4");
	
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
	std::cout << "DONE (" << timestepCount << ")\n";
	
	size_t passCount = (_solutionInterval==0) ? 1 : (timestepCount + _solutionInterval - 1) / _solutionInterval;
	
	if(!_modelFilename.empty()) {
			std::cout << "Reading model... " << std::flush;
			_model = Model(_modelFilename.c_str());
			std::cout << "DONE\n";
	}

	if(_saveSolutionsFiles)
	{
		for(size_t ant=0; ant!=antennaCount; ++ant)
		{
			std::ostringstream antFilename;
			antFilename << "peel-sol-ant" << ant << ".txt";
			std::ofstream(antFilename.str().c_str());
		}
	}
		
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
			calMethods[ch]->SetOnlySolveScalar(_onlyScalar);
			calMethods[ch]->SetOnlySolveDiag(_onlyDiag);
			calMethods[ch]->SetOnlySolveRotation(_onlyRotation);
		}
		std::unique_ptr<MSPredicter> predicter;
		std::unique_ptr<BeamEvaluator> beamEvaluator;
		std::vector<std::complex<double>> beamValues;
		if(_model.Empty()) {
			std::cout << "Reading data and model column... " << std::flush;
			modelColumn.reset(new casacore::ROArrayColumn<complex_t>(_ms, _ms.columnName(casacore::MSMainEnums::MODEL_DATA)));
		}
		else {
			if(_beamOnSource || _applyBeam)
			{
				beamEvaluator.reset(new BeamEvaluator(_ms, true, ""));
			}
			if(_beamOnSource)
			{
				if(_model.SourceCount() != 1)
					std::cout << "Warning: To correct for the beam, there should be exactly one source in the model";
				const ModelComponent& component = _model.Source(0).Peak();
				beamValues.resize(channelCount*4);
				double beamSum[4] = {0.0, 0.0, 0.0, 0.0};
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					double frequency = _bandData.ChannelFrequency(ch);
					beamEvaluator->EvaluateApparentToAbsGain(component.PosRA(), component.PosDec(), frequency, &beamValues[ch*4]);
					for(size_t p=0; p!=4; ++p)
						beamSum[p] += std::abs(beamValues[ch*4+p]);
				}
			}
			
			predicter.reset(new MSPredicter(_ms, _threadCount, _model));
			predicter->SetApplyBeam(_applyBeam);
			predicter->SetStartRow(startRow);
			predicter->SetEndRow(endRow);
		}
		
		std::vector<std::complex<double> > modelValues(4 * channelCount);
		casacore::Array<complex_t> data(dataShape), modelData(dataShape);
		casacore::Array<float> weights(dataShape);
		casacore::Array<bool> flags(dataShape);
		time = timeColumn(startRow);
		size_t selectedCount = 0, notSelected = 0;
		MSPredicter::RowData rowData;
		
		predicter->Start(false);
		while(predicter->GetNextRow(rowData))
		{
			// Cross correlation?
			size_t antenna1 = rowData.a1, antenna2 = rowData.a2;
			if(antenna1 == antenna2)
			{
			}
			else {
				size_t rowIndex = rowData.rowIndex;
				std::unique_lock<std::mutex> lock(predicter->IOMutex());
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
					size_t timeIndex = rowData.timeIndex;
					for(size_t ch = 0; ch!=channelCount; ++ch)
					{
						size_t chIndex = ch * 4;
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
			threadGroup.add_thread(new boost::thread(&Peeler::calibrateThreadFunction, this, threadData));
		}
		threadGroup.join_all();

		if(_saveSolutionsFiles)
		{
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
						val[p] = calMethods[ch]->JonesSolution(ant, p);
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
		}
			
		/**
		 * Do the subtraction
		 */
		predicter.reset(new MSPredicter(_ms, _threadCount, _model));
		predicter->SetApplyBeam(_applyBeam);
		predicter->SetStartRow(startRow);
		predicter->SetEndRow(endRow);
		predicter->Start(false);
		while(predicter->GetNextRow(rowData))
		{
			size_t rowIndex = rowData.rowIndex;
			size_t antenna1 = rowData.a1, antenna2 = rowData.a2;
			if(antenna1 != antenna2)
			{
				std::unique_lock<std::mutex> lock(predicter->IOMutex());
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				lock.unlock();
				casacore::Complex *dataPtr = data.cbegin();
				bool *flagPtr = flags.cbegin();
				
				// Apply solutions to model and subtract from data
				bool flagsWereChanged = false;
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					std::complex<double> solutionsA[4], solutionsB[4];
					for(size_t p=0; p!=4; ++p)
					{
						solutionsA[p] = calMethods[ch]->JonesSolution(antenna1, p);
						solutionsB[p] = calMethods[ch]->JonesSolution(antenna2, p);
					}
					if(Matrix2x2::IsFinite(solutionsA) && Matrix2x2::IsFinite(solutionsB))
					{
						std::complex<double> temp[4], doubleData[4];
						Matrix2x2::ATimesB(temp, solutionsA, rowData.modelData+ch*4);
						// Store temporarily in solutionsA
						Matrix2x2::ATimesHermB(solutionsA, temp, solutionsB);
						Matrix2x2::Assign(doubleData, dataPtr+ch*4);
						Matrix2x2::Subtract(doubleData, solutionsA);
						Matrix2x2::Assign(dataPtr+ch*4, doubleData);
					}
					else {
						for(size_t p=0; p!=4; ++p)
							flagPtr[ch*4 + p] = true;
						flagsWereChanged = true;
					}
				}
				
				lock.lock();
				dataColumn.put(rowIndex, data);
				if(flagsWereChanged)
					flagColumn.put(rowIndex, flags);
			}
			
			predicter->FinishRow(rowData);
		}
		
		/**
		 * Deallocate resources
		 */
		for(size_t ch=0; ch!=channelCount; ++ch)
			delete calMethods[ch];
	}
}

void Peeler::calibrateThreadFunction(Peeler::ThreadData data)
{
	boost::mutex::scoped_lock lock(*data.mutex);
	size_t lastSuccessfulChannel = data.tasks->front();
	while(!data.tasks->empty()) {
		size_t taskIndex = data.tasks->front();
		data.tasks->pop();
		lock.unlock();
		if(lastSuccessfulChannel != taskIndex)
			(*(data.calMethods))[taskIndex]->InitSolutions(*(*(data.calMethods))[lastSuccessfulChannel]);
		size_t iters = _nIter;
		double limit = _stoppingAccuracy;
		(*(data.calMethods))[taskIndex]->Execute(limit, iters);
		if((iters >= _nIter || !std::isfinite(limit)) && !(*(data.calMethods))[taskIndex]->OnlySolveRotation())
		{
			std::cout << "Recalculating channel " << taskIndex << " (accuracy=" << limit << ").\n";
			(*(data.calMethods))[taskIndex]->InitSolutionsToUnity();
			iters = _nIter;
			limit = _stoppingAccuracy;
			(*(data.calMethods))[taskIndex]->Execute(limit, iters);

			if((iters >= _nIter && limit > _minAccuracy) || !std::isfinite(limit))
			{
				std::cout << "Channel " << taskIndex << " did not converge (accuracy=" << limit << "), setting gains to NaN.\n";
				(*(data.calMethods))[taskIndex]->InitSolutionsToNaN();
			}
			else {
				if(iters >= _nIter && limit > _stoppingAccuracy)
				{
					std::cout << "Channel " << taskIndex << " converged (accuracy=" << limit << ") but did not reach stopping accuracy.\n";
				}
				lastSuccessfulChannel = taskIndex;
			}
		}
		else {
			lastSuccessfulChannel = taskIndex;
		}
		lock.lock();
	}
}
