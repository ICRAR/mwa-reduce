#include "ionpeeler.h"

#include "beamevaluator.h"
#include "ionpeeler.h"
#include "imageweights.h"
#include "predicter.h"
#include "progressbar.h"
#include "serializable.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#ifdef HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#endif

IonPeeler::IonPeeler() : _solutionInterval(1), _fitIterationCount(3), _applyBeam(true), _weightMode(WeightMode::NaturalWeighted), _weightGridSize(0), _weightPixelScale(0.0)
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
	_model.SortOnBrightness();
	
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
	_antennaCount = ms.antenna().nrow();
	
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
		
		_predicters[s] = new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), channelCount);
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
	_passCount = (_solutionInterval==0) ? 1 : (timestepCount + _solutionInterval - 1) / _solutionInterval;
	casa::Array<std::complex<float> > data(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	_solutions.resize(_model.SourceCount() * _passCount * channelCount);
	_solutionWeights.resize(_model.SourceCount() * _passCount * channelCount);
	_dataArrays.resize(channelCount);
	_weightArrays.resize(channelCount);
	for(_pass=0; _pass!=_passCount; ++_pass)
	{
		_startTimestep = timestepCount * _pass / _passCount;
		_endTimestep = timestepCount * (_pass+1) / _passCount;
		const size_t timestepsInPass = _endTimestep - _startTimestep;
		_curStartRow = timestepRows[_startTimestep];
		_curEndRow = timestepRows[_endTimestep];
		_rowData.resize(_curEndRow - _curStartRow);
		
		if(beamEvaluator.Time().getValue() != timeMColumn(_curStartRow).getValue())
		{
			std::cout << "Evaluating beam...\n";
			beamEvaluator.SetTime(timeMColumn(_curStartRow));
			for(size_t s=0; s!=_model.SourceCount(); ++s)
				_predicters[s]->UpdateBeam(_predictionModels[s]);
		}
		
		std::cout << "Reading (T " << _startTimestep << "-" << _endTimestep << ", rows " << _curStartRow << '-' << _curEndRow << ")...\n";
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			_dataArrays[ch] = new VisibilityArray<std::complex<double>, 4>(1, _antennaCount, timestepsInPass);
			_weightArrays[ch] = new VisibilityArray<double, 2>(1, _antennaCount, timestepsInPass);
		}
		double time = timeColumn(_curStartRow);
		size_t timeIndex = _startTimestep;
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
					
					std::complex<double> *arrPtr = _dataArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep);
					arrPtr[0] = dataPtr[0]; arrPtr[1] = dataPtr[1];
					arrPtr[2] = dataPtr[2]; arrPtr[3] = dataPtr[3];
					double* weightsWritePtr = _weightArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep);
					weightsWritePtr[0] = weightsReadPtr[0];
					weightsWritePtr[1] = weightsReadPtr[3];
					
					weightsReadPtr += 4;
					dataPtr += 4;
					flagPtr += 4;
				}
			}
		}
		
		_progressBar.reset(new ProgressBar("Processing channels"));
		_stats = PeelingStats();
		std::vector<size_t> tasks;
		for(size_t ch=0; ch!=channelCount; ++ch)
			tasks.push_back(channelCount - ch - 1);
		std::mutex mutex;
		std::vector<std::thread> threads;
		for(size_t i=0; i!=_cpuCount; ++i)
			threads.push_back(std::thread(&IonPeeler::processingThreadFunction, this, &mutex, &tasks));
		for(std::vector<std::thread>::iterator t = threads.begin(); t!=threads.end(); ++t)
			t->join();
		
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
					std::complex<double> *arrPtr = _dataArrays[ch]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep);
					dataPtr[0] = arrPtr[0]; dataPtr[1] = arrPtr[1];
					dataPtr[2] = arrPtr[2]; dataPtr[3] = arrPtr[3];
					dataPtr += 4;
				}
				dataColumn.put(rowIndex, data);
			}
		}
		_progressBar.reset();
		
		outputStats(_stats);
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			delete _dataArrays[ch];
			delete _weightArrays[ch];
		}
	}
	
	for(size_t s=0; s!=_model.SourceCount(); ++s)
		_predicters[s] = new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), channelCount);
}

void IonPeeler::processingThreadFunction(std::mutex* mutex, std::vector<size_t>* tasks)
{
	PeelingStats stats;
	std::unique_lock<std::mutex> lock(*mutex);
	while(!tasks->empty())
	{
		_progressBar->SetProgress(_bandData.ChannelCount() - tasks->size(), _bandData.ChannelCount());
		const size_t channel = tasks->back();
		tasks->pop_back();
		lock.unlock();
		
		processChannel(channel, stats);
		
		lock.lock();
	}
	// (mutex is still locked)
	_stats += stats;
}

void IonPeeler::processChannel(size_t channelIndex, PeelingStats& stats)
{
	//scalarGainFitter(channelIndex);
	positionFitter(channelIndex, stats);
}

int IonPeeler::posMinimizationFunc(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	const FittingInfo &fittingInfo = *reinterpret_cast<FittingInfo*>(data);
	const IonPeeler &ionPeeler = *reinterpret_cast<IonPeeler*>(fittingInfo.ionPeeler);
	const size_t startTimestep = ionPeeler._startTimestep;
	const double lambda = ionPeeler._bandData.ChannelWavelength(fittingInfo.channelIndex);
	double g = gsl_vector_get(xvec, 0);
	double dl = gsl_vector_get(xvec, 1);
	double dm = gsl_vector_get(xvec, 2);
	
	size_t dataPointIndex = 0;
	const std::complex<double>* modelPtr = &fittingInfo.modelData[0];
	for(size_t row=ionPeeler._curStartRow; row!=ionPeeler._curEndRow; ++row)
	{
		const RowData& rowData = ionPeeler._rowData[row - ionPeeler._curStartRow];
		if(rowData.a1 != rowData.a2)
		{
			const std::complex<double>* dataPtr = ionPeeler._dataArrays[fittingInfo.channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimestep);
			std::complex<double> value = 0.0;
			if(isfinite(dataPtr[0]) && isfinite(dataPtr[3]))
			{
				double uOverL = rowData.u / lambda, vOverL = rowData.v / lambda;
				double expTerm = (uOverL * dl + vOverL * dm) * 2.0 * M_PI;
				double sExp, cExp;
				sincos(expTerm, &sExp, &cExp);
				std::complex<double> rotModelXX = std::complex<double>(
						modelPtr[0].real() * cExp - modelPtr[0].imag() * sExp,
						modelPtr[0].real() * sExp + modelPtr[0].imag() * cExp
					);
				std::complex<double> rotModelYY = std::complex<double>(
						modelPtr[3].real() * cExp - modelPtr[3].imag() * sExp,
						modelPtr[3].real() * sExp + modelPtr[3].imag() * cExp
					);
				value = (g * (rotModelXX + rotModelYY) - dataPtr[0] - dataPtr[3]) * 0.5;
			}
			gsl_vector_set(f, dataPointIndex, value.real());
			++dataPointIndex;
			gsl_vector_set(f, dataPointIndex, value.imag());
			++dataPointIndex;
			modelPtr += 4;
		}
	}
		
	return GSL_SUCCESS;
}
int IonPeeler::posMinimizationFuncDeriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	const FittingInfo &fittingInfo = *reinterpret_cast<FittingInfo*>(data);
	const IonPeeler &ionPeeler = *reinterpret_cast<IonPeeler*>(fittingInfo.ionPeeler);
	const size_t startTimestep = ionPeeler._startTimestep;
	const double lambda = ionPeeler._bandData.ChannelWavelength(fittingInfo.channelIndex);
	
	double g = gsl_vector_get(xvec, 0);
	double dl = gsl_vector_get(xvec, 1);
	double dm = gsl_vector_get(xvec, 2);
	
	size_t dataPointIndex = 0;
	const std::complex<double>* modelPtr = &fittingInfo.modelData[0];
	for(size_t row=ionPeeler._curStartRow; row!=ionPeeler._curEndRow; ++row)
	{
		const RowData& rowData = ionPeeler._rowData[row - ionPeeler._curStartRow];
		if(rowData.a1 != rowData.a2)
		{
			std::complex<double> dfdg = 0.0, dfddl = 0.0, dfddm = 0.0;
			const std::complex<double>* dataPtr = ionPeeler._dataArrays[fittingInfo.channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimestep);
			if(isfinite(dataPtr[0]) && isfinite(dataPtr[3]))
			{
				double uOverL = rowData.u / lambda, vOverL = rowData.v / lambda;
				double expTerm = (uOverL * dl + vOverL * dm) * 2.0 * M_PI;
				double sExp, cExp;
				sincos(expTerm, &sExp, &cExp);
				std::complex<double> rotModelXX = std::complex<double>(
						modelPtr[0].real() * cExp - modelPtr[0].imag() * sExp,
						modelPtr[0].real() * sExp + modelPtr[0].imag() * cExp
					);
				std::complex<double> rotModelYY = std::complex<double>(
						modelPtr[3].real() * cExp - modelPtr[3].imag() * sExp,
						modelPtr[3].real() * sExp + modelPtr[3].imag() * cExp
					);
				dfdg = (rotModelXX + rotModelYY) * 0.5;
				std::complex<double> sum = rotModelXX + rotModelYY;
				sum = std::complex<double>(-sum.imag(), sum.real());
				dfddl = sum * g * rowData.u * M_PI;
				dfddm = sum * g * rowData.v * M_PI;
			}
			gsl_matrix_set(J, dataPointIndex, 0, dfdg.real());
			gsl_matrix_set(J, dataPointIndex, 1, dfddl.real());
			gsl_matrix_set(J, dataPointIndex, 2, dfddm.real());
			++dataPointIndex;
			gsl_matrix_set(J, dataPointIndex, 0, dfdg.imag());
			gsl_matrix_set(J, dataPointIndex, 1, dfddl.imag());
			gsl_matrix_set(J, dataPointIndex, 2, dfddm.imag());
			++dataPointIndex;
			modelPtr += 4;
		}
	}
	return GSL_SUCCESS;
}
int IonPeeler::posMinimizationFuncBoth(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	posMinimizationFunc(x, data, f);
	posMinimizationFuncDeriv(x, data, J);
	return GSL_SUCCESS;
}

void IonPeeler::positionFitter(size_t channelIndex, PeelingStats& stats)
{
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	const double lambda = _bandData.ChannelWavelength(channelIndex);
	size_t dataValueCount = 0;
	for(size_t row=_curStartRow; row!=_curEndRow; ++row) 
	{
		const RowData& rowData = _rowData[row - _curStartRow];
		if(rowData.a1 != rowData.a2) ++dataValueCount;
	}
	dataValueCount *= 2; // For real and imaginary
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(T, dataValueCount, 3);
	
	FittingInfo fInfo;
	fInfo.ionPeeler = this;
	fInfo.modelData.resize((_curEndRow-_curStartRow)*4);
	fInfo.channelIndex = channelIndex;
	
	ao::uvector<double>
		_solutionsG(_predictionModels.size()),
		_solutionsDL(_predictionModels.size()),
		_solutionsDM(_predictionModels.size());
	
	for(size_t fitIteration=0; fitIteration!=_fitIterationCount; ++fitIteration)
	{
		for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
		{
			double&
				g = _solutionsG[sourceIndex],
				dl = _solutionsDL[sourceIndex],
				dm = _solutionsDM[sourceIndex];
			
			// Predict
			ao::uvector<std::complex<double>>::iterator modelPtr = fInfo.modelData.begin();
			for(size_t row=_curStartRow; row!=_curEndRow; ++row)
			{
				const RowData& rowData = _rowData[row - _curStartRow];
				if(rowData.a1 != rowData.a2)
				{
					const double uOverL = rowData.u/lambda, vOverL = rowData.v/lambda, wOverL = rowData.w/lambda;
					_predicters[sourceIndex]->Predict4(&*modelPtr, _predictionModels[sourceIndex], uOverL, vOverL, wOverL, channelIndex, rowData.a1, rowData.a2);
					modelPtr += 4;
				}
			}
			
			// Add back if it was subtracted
			if(fitIteration != 0)
			{
				modelPtr = fInfo.modelData.begin();
				for(size_t row=_curStartRow; row!=_curEndRow; ++row)
				{
					const RowData& rowData = _rowData[row - _curStartRow];
					if(rowData.a1 != rowData.a2)
					{
						std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - _startTimestep);
						double uOverL = rowData.u / lambda, vOverL = rowData.v / lambda;
						double expTerm = (uOverL * dl + vOverL * dm) * 2.0 * M_PI;
						double sExp, cExp;
						sincos(expTerm, &sExp, &cExp);
						for(size_t p=0; p!=4; ++p)
						{
							std::complex<double> rotModel = std::complex<double>(
									modelPtr[p].real() * cExp - modelPtr[p].imag() * sExp,
									modelPtr[p].real() * sExp + modelPtr[p].imag() * cExp
								);
							std::complex<double> corModel = g * rotModel;
							dataPtr[p] -= corModel;
						}
						modelPtr += 4;
					}
				}
			}
			
			// Fit
			gsl_multifit_function_fdf fdf;
			fdf.f = &posMinimizationFunc;
			fdf.df = &posMinimizationFuncDeriv;
			fdf.fdf = &posMinimizationFuncBoth;
			fdf.n = dataValueCount;
			fdf.p = 3;
			fdf.params = &fInfo;
			
			double initialValsArray[3];
			if(fitIteration == 0)
			{
				initialValsArray[0] = 1.0;
				initialValsArray[1] = 0.0;
				initialValsArray[2] = 0.0;
			}
			else {
				initialValsArray[0] = g;
				initialValsArray[1] = dl;
				initialValsArray[2] = dm;
			}
			gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 3);
			gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

			int status;
			size_t iter = 0;
			do {
				iter++;
				status = gsl_multifit_fdfsolver_iterate(solver);
				
				if(status)
					break;
				
				status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
				
			} while (status == GSL_CONTINUE && iter < 100);
			
			g = gsl_vector_get (solver->x, 0);
			dl = gsl_vector_get (solver->x, 1);
			dm = gsl_vector_get (solver->x, 2);
			//std::cout << gsl_strerror(status) << "\n";
			//std::cout << (status==GSL_CONTINUE ? "CONTINUE " : "X ") << iter << ", g=" << g << ",dl=" << dl << ",dm=-" << dm << '\n';
			
			if(fitIteration+1 == _fitIterationCount)
			{
				stats.lsFits++;
				stats.lsFittingIterations += iter;
				stats.gSum += g;
				stats.dlSum += dl;
				stats.dmSum += dm;
				stats.gSumSq += g*g;
				stats.dlSumSq += dl*dl;
				stats.dmSumSq += dm*dm;
			}
			
			// Subtract
			modelPtr = fInfo.modelData.begin();
			for(size_t row=_curStartRow; row!=_curEndRow; ++row)
			{
				const RowData& rowData = _rowData[row - _curStartRow];
				if(rowData.a1 != rowData.a2)
				{
					std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - _startTimestep);
					double uOverL = rowData.u / lambda, vOverL = rowData.v / lambda;
					double expTerm = (uOverL * dl + vOverL * dm) * 2.0 * M_PI;
					double sExp, cExp;
					sincos(expTerm, &sExp, &cExp);
					for(size_t p=0; p!=4; ++p)
					{
						std::complex<double> rotModel = std::complex<double>(
								modelPtr[p].real() * cExp - modelPtr[p].imag() * sExp,
								modelPtr[p].real() * sExp + modelPtr[p].imag() * cExp
							);
						std::complex<double> corModel = g * rotModel;
						dataPtr[p] -= corModel;
					}
					modelPtr += 4;
				}
			}
		}
	}
	gsl_multifit_fdfsolver_free(solver);
}

void IonPeeler::scalarGainFitter(size_t channelIndex)
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
					_predicters[sourceIndex]->Predict4(&*modelPtr, _predictionModels[sourceIndex], uOverL, vOverL, wOverL, channelIndex, rowData.a1, rowData.a2);
					
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
			_solutions[(_pass * _bandData.ChannelCount() + channelIndex)*_predictionModels.size() + sourceIndex] = ionTerm;
			_solutionWeights[(_pass * _bandData.ChannelCount() + channelIndex)*_predictionModels.size() + sourceIndex] = ionWeight;
			
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
					Matrix2x2::ScalarMultiply(&*modelPtr, ionTerm);
					
					std::complex<double>* dataPtr = _dataArrays[channelIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimeIndex);
					Matrix2x2::Subtract(dataPtr, &*modelPtr);
					modelPtr += 4;
				}
			}
		}
		else {
			_solutions[(_pass * _bandData.ChannelCount() + channelIndex)*_predictionModels.size() + sourceIndex] = 0.0;
			_solutionWeights[(_pass * _bandData.ChannelCount() + channelIndex)*_predictionModels.size() + sourceIndex] = 0.0;
		}
	}
}


void IonPeeler::SaveSolutions(const string& filename) const
{
	std::ofstream str(filename.c_str());
	Serializable::SerializeToUInt64(str, _passCount);
	Serializable::SerializeToUInt64(str, _bandData.ChannelCount());
	Serializable::SerializeToUInt64(str, _predictionModels.size());
	str.write(reinterpret_cast<const char*>(_solutions.data()), sizeof(std::complex<double>) * _solutions.size());
	str.write(reinterpret_cast<const char*>(_solutionWeights.data()), sizeof(double) * _solutionWeights.size());
}

std::string IonPeeler::radToString(double r)
{
	std::ostringstream str;
	str << 60.0*60.0*(r*180.0)/M_PI << "''";
	return str.str();
}

void IonPeeler::outputStats(const IonPeeler::PeelingStats& stats)
{
	std::cout
		<< "Avg lsiterations=" << double(stats.lsFittingIterations) / stats.lsFits
		<< ", avg g=" << stats.gSum / stats.lsFits
		<< ", ion avg=" << radToString(stats.dlSum/stats.lsFits) << ',' << radToString(stats.dmSum/stats.lsFits)
		<< ", ion RMS=" << radToString(sqrt(stats.dlSumSq/stats.lsFits)) << ',' << radToString(sqrt(stats.dmSumSq/stats.lsFits))
		<< '\n';
}
