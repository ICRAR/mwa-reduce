#include "ionpeeler.h"

#include "beamevaluator.h"
#include "ionpeeler.h"
#include "imageweights.h"
#include "predicter.h"
#include "progressbar.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

IonPeeler::IonPeeler() : _solutionInterval(1), _applyBeam(true), _weightMode(WeightMode::NaturalWeighted), _weightGridSize(0), _weightPixelScale(0.0)
{ }

IonPeeler::~IonPeeler()
{ }

void IonPeeler::initWeighting(casa::MeasurementSet& ms)
{
	if(_weightMode.RequiresGridding())
	{
		std::cout << "Precalculating weights for " << _weightMode.ToString() << " weighting...\n";
		_imageWeights.reset(new ImageWeights(_weightGridSize, _weightGridSize, _weightPixelScale));
		_imageWeights->Grid(ms, _weightMode);
	}
}

void IonPeeler::Peel(const char* msName, const char* modelName)
{
	casa::MeasurementSet ms(msName, casa::MeasurementSet::Update);
	_model = Model(modelName);
	//_model.SortOnBrightness();
	
	initWeighting(ms);
	
	size_t startRow = 0;
	std::string solutionFile;
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ArrayColumn<casa::Complex> dataColumn(ms, _dataColumnName);
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	_bandData = BandData(ms.spectralWindow());
	const size_t channelCount = _bandData.ChannelCount();
	const size_t antennaCount = ms.antenna().nrow();
	
	casa::MSField fieldTable = ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	// By setting the time beforehand, we don't waste time calculating a time step we don't need.
	casa::MEpoch::ROScalarColumn timeMColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	casa::MEpoch startTime = timeMColumn(startRow);
	BeamEvaluator beamEvaluator(ms);
	beamEvaluator.SetTime(startTime);
	
	std::cout << "Initializing model predicter for " << _model.SourceCount() << " sources...\n";
	_predictionModels.resize(_model.SourceCount());
	_predicters.resize(_model.SourceCount());
	for(size_t s=0; s!=_model.SourceCount(); ++s)
	{
		_predictionModels[s].AddSource(_model.Source(s));
		
		_predicters[s].reset(new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), channelCount));
		if(_applyBeam)
			_predicters[s]->Initialize(_predictionModels[s], solutionFile, &beamEvaluator);
		else
			_predicters[s]->Initialize(_predictionModels[s], solutionFile);
	}
	_predicters.front()->ReportSources(_predictionModels.front());
	
	std::cout << "Counting timesteps...\n";
	double time = -1.0;
	ao::uvector<size_t> timestepRows;
	casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
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
	
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	
	_cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	size_t passCount = (_solutionInterval==0) ? 1 : (timestepCount + _solutionInterval - 1) / _solutionInterval;
	casa::Array<std::complex<float> > data(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	_dataArrays.resize(channelCount);
	_weightArrays.resize(channelCount);
	for(size_t pass=0; pass!=passCount; ++pass)
	{
		const size_t
			startTimestep = timestepCount * pass / passCount,
			endTimestep = timestepCount * (pass+1) / passCount,
			timestepsInPass = endTimestep - startTimestep;
		_curStartRow = timestepRows[startTimestep];
		_curEndRow = timestepRows[endTimestep];
		_rowData.resize(_curEndRow - _curStartRow);
		
		if(beamEvaluator.Time().getValue() != timeMColumn(_curStartRow).getValue())
		{
			std::cout << "Evaluating beam...\n";
			beamEvaluator.SetTime(timeMColumn(_curStartRow));
			for(size_t s=0; s!=_model.SourceCount(); ++s)
				_predicters[s]->UpdateBeam(_predictionModels[s]);
		}
		
		std::cout << "Reading (T " << startTimestep << "-" << endTimestep << ", rows " << _curStartRow << '-' << _curEndRow << ")...\n";
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			_dataArrays[ch].reset(new VisibilityArray<std::complex<double>, 4>(1, antennaCount, timestepsInPass));
			_weightArrays[ch].reset(new VisibilityArray<double, 2>(1, antennaCount, timestepsInPass));
		}
		double time = timeColumn(_curStartRow);
		size_t timeIndex = startTimestep;
		for(size_t rowIndex=_curStartRow; rowIndex!=_curEndRow; ++rowIndex)
		{
			RowData& rowData = _rowData[rowIndex-_curStartRow];
			rowData.a1 = ant1Column(rowIndex);
			rowData.a2 = ant2Column(rowIndex);
			if(timeColumn(rowIndex) != time)
			{
				++timeIndex;
				time = timeColumn(rowIndex);
			}
			rowData.timeIndex = timeIndex;
			if(rowData.a1 != rowData.a2)
			{
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
				rowData.u = *uvwI; ++uvwI;
				rowData.v = *uvwI; ++uvwI;
				rowData.w = *uvwI;
			
				dataColumn.get(rowIndex, data);
				weightColumn.get(rowIndex, weights);
				flagColumn.get(rowIndex, flags);
				std::complex<float> *dataPtr = data.cbegin();
				float *weightsReadPtr = weights.cbegin();
				bool *flagPtr = flags.cbegin();
				
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					for(size_t p=0; p!=4; ++p)
					{
						if(flagPtr[p]) weightsReadPtr[p] = 0.0;
					}
					
					std::complex<double> *arrPtr = _dataArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - startTimestep);
					arrPtr[0] = dataPtr[0]; arrPtr[1] = dataPtr[1];
					arrPtr[2] = dataPtr[2]; arrPtr[3] = dataPtr[3];
					double* weightsWritePtr = _weightArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - startTimestep);
					weightsWritePtr[0] = weightsReadPtr[0];
					weightsWritePtr[1] = weightsReadPtr[3];
					
					weightsReadPtr += 4;
					dataPtr += 4;
					flagPtr += 4;
				}
			}
		}
		
		_progressBar.reset(new ProgressBar("Processing channels"));
		std::vector<size_t> tasks;
		for(size_t ch=0; ch!=channelCount; ++ch)
			tasks.push_back(channelCount - ch - 1);
		std::mutex mutex;
		std::vector<std::thread> threads;
		for(size_t i=0; i!=_cpuCount; ++i)
			threads.push_back(std::thread(&IonPeeler::processingThreadFunction, this, &mutex, &tasks));
		for(auto& t : threads)
			t.join();
		
		// Write back
		for(size_t rowIndex=_curStartRow; rowIndex!=_curEndRow; ++rowIndex)
		{
			RowData& rowData = _rowData[rowIndex-_curStartRow];
			if(rowData.a1 != rowData.a2)
			{
				timeIndex = rowData.timeIndex;
				std::complex<float> *dataPtr = data.cbegin();
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					std::complex<double> *arrPtr = _dataArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - startTimestep);
					dataPtr[0] = arrPtr[0]; dataPtr[1] = arrPtr[1];
					dataPtr[2] = arrPtr[2]; dataPtr[3] = arrPtr[3];
					dataPtr += 4;
				}
				dataColumn.put(rowIndex, data);
			}
		}
		_progressBar.reset();
	}
}

void IonPeeler::processingThreadFunction(std::mutex* mutex, std::vector<size_t>* tasks)
{
	std::unique_lock<std::mutex> lock(*mutex);
	while(!tasks->empty())
	{
		_progressBar->SetProgress(_bandData.ChannelCount() - tasks->size(), _bandData.ChannelCount());
		const size_t channel = tasks->back();
		tasks->pop_back();
		lock.unlock();
		
		processChannel(channel);
		
		lock.lock();
	}
}

void IonPeeler::processChannel(size_t channelIndex)
{
	double lambda = _bandData.ChannelWavelength(channelIndex);
	size_t startTimeIndex = _rowData.front().timeIndex;
	
	// Pre-measure the weight: this is equal for all sources so do only once
	double ionWeight = 0.0;
	for(size_t row=_curStartRow; row!=_curEndRow; ++row)
	{
		const RowData& rowData = _rowData[row - _curStartRow];
		if(rowData.a1 != rowData.a2)
		{
			const std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
			const double* weightPtr = _weightArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
			double weight;
			if(isfinite(dataPtr[0]) && isfinite(dataPtr[3]) && dataPtr[0] != 0.0 && dataPtr[3] != 0.0)
				weight = weightPtr[0] + weightPtr[1];
			else
				weight = 0.0;
			
			if(_imageWeights != 0)
			{
				const double uOverL = rowData.u/lambda, vOverL = rowData.v/lambda;
				switch(_weightMode.Mode())
				{
					case WeightMode::UniformWeighted:
						weight *= _imageWeights->GetUniformWeight(uOverL, vOverL);
						break;
					case WeightMode::BriggsWeighted:
						weight *= _imageWeights->GetBriggsWeight(uOverL, vOverL);
						break;
					case WeightMode::NaturalWeighted:
					case WeightMode::DistanceWeighted:
						break;
				}
			}
			ionWeight += weight;
		}
	}
	
	// Measure ionospheric term
	for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
	{
		std::complex<double> ionTerm = 0.0;
		ao::uvector<std::complex<double>> modelData((_curEndRow-_curStartRow)*4);
		ao::uvector<std::complex<double>>::iterator modelPtr = modelData.begin();
		for(size_t row=_curStartRow; row!=_curEndRow; ++row)
		{
			const RowData& rowData = _rowData[row - _curStartRow];
			if(rowData.a1 != rowData.a2)
			{
				const double uOverL = rowData.u/lambda, vOverL = rowData.v/lambda, wOverL = rowData.w/lambda;
				const std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
				const double* weightPtr = _weightArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
				
				if(isfinite(dataPtr[0]) && isfinite(dataPtr[3]) && dataPtr[0] != 0.0 && dataPtr[3] != 0.0)
				{
					// Predict visibility
					_predicters[sourceIndex]->Predict4(modelPtr, _predictionModels[sourceIndex], uOverL, vOverL, wOverL, channelIndex, rowData.a1, rowData.a2);
					
					double imageWeight = 1.0;
					if(_imageWeights != 0)
					{
						switch(_weightMode.Mode())
						{
							case WeightMode::UniformWeighted:
								imageWeight = _imageWeights->GetUniformWeight(uOverL, vOverL);
								break;
							case WeightMode::BriggsWeighted:
								imageWeight = _imageWeights->GetBriggsWeight(uOverL, vOverL);
								break;
							case WeightMode::NaturalWeighted:
							case WeightMode::DistanceWeighted:
								break;
						}
					}
					
					// Calculate ionospheric term
					ionTerm += (weightPtr[0] * dataPtr[3] / modelPtr[0] + weightPtr[1] * dataPtr[3] / modelPtr[3]) * imageWeight;
				}
				modelPtr += 4;
			}
		}
		
		if(ionWeight != 0.0 && isfinite(ionTerm))
		{
			ionTerm /= ionWeight;
			
			/*if((channelIndex + 11) % 16 == 0)
			{
				std::cout << "T=" << startTimeIndex << ", ch=" << channelIndex << ", gain=" << std::abs(ionTerm) << ", phase=" << std::atan2(ionTerm.imag(), ionTerm.real()) << " (w=" << ionWeight << ")\n";
			}*/
		
			// Subtract source from visibilities in mem
			modelPtr = modelData.begin();
			for(size_t row=_curStartRow; row!=_curEndRow; ++row)
			{
				const RowData& rowData = _rowData[row - _curStartRow];
				if(rowData.a1 != rowData.a2)
				{
					Matrix2x2::ScalarMultiply(modelPtr, ionTerm);
					
					std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
					Matrix2x2::Subtract(dataPtr, modelPtr);
					modelPtr += 4;
				}
			}
		}
	}
}
