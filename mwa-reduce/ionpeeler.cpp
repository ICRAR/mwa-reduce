#include "ionpeeler.h"

#include "beamevaluator.h"
#include "ionpeeler.h"
#include "ionsolutionfile.h"
#include "imageweights.h"
#include "predicter.h"
#include "progressbar.h"
#include "serializable.h"
#include "msselection.h"
#include "imagecoordinates.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <boost/thread/thread.hpp>

#ifdef HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#endif

IonPeeler::IonPeeler() :
	_solutionInterval(1), _fitIterationCount(3), _applyBeam(true),
	_channelBlockSize(0), _channelBlockCount(1),
	_weightMode(WeightMode::NaturalWeighted), _weightGridSize(0), _weightPixelScale(0.0),
	_clusterFluxLimit(0.0), _distanceLimit(0.0), _verbose(false)
{ }

IonPeeler::~IonPeeler()
{ }

void IonPeeler::initWeighting(casa::MeasurementSet& ms)
{
	if(_weightMode.RequiresGridding())
	{
		std::cout << "Precalculating weights for " << _weightMode.ToString() << " weighting...\n";
		_imageWeights.reset(new ImageWeights(_weightGridSize, _weightGridSize, _weightPixelScale, _weightPixelScale));
		_imageWeights->Grid(ms, _weightMode, MSSelection::Everything());
	}
}

void IonPeeler::Peel(const char* msName, const char* modelName, const char* solutionFilename)
{
	casa::MeasurementSet ms(msName, casa::MeasurementSet::Update);
	
	if(_dataColumnName.empty())
	{
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			std::cout << "Measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_dataColumnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in measurement set: tasks will be applied on the data column.\n";
			_dataColumnName= "DATA";
		}
	}
	
	_model = Model(modelName);
	
	initWeighting(ms);
	
	size_t startRow = 0;
	
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
	
	std::vector<std::string> clusterNames;
	_model.GetClusterNames(clusterNames);
	std::cout << "Requesting to solve for " << clusterNames.size() << " directions.\n";
	_predictionModels.resize(clusterNames.size());
	_predicters.resize(clusterNames.size());
	
	std::set<size_t> clustersToRemove;
	for(size_t clusterIter=0; clusterIter!=clusterNames.size(); ++clusterIter)
	{
		SourceGroup cluster;
		_model.GetSourcesInCluster(clusterNames[clusterIter], cluster);
		for(std::vector<ModelSource>::const_iterator sourceIter=cluster.begin();
				sourceIter!=cluster.end(); ++sourceIter)
		{
			_predictionModels[clusterIter].AddSource(*sourceIter);
		}
		
		_predicters[clusterIter] = new Predicter(phaseCentreRA, phaseCentreDec, _bandData.LowestFrequency(), _bandData.HighestFrequency(), channelCount);
		if(_applyBeam)
			_predicters[clusterIter]->Initialize(_predictionModels[clusterIter], "", &beamEvaluator);
		else
			_predicters[clusterIter]->Initialize(_predictionModels[clusterIter]);
		
		double distance =
			ImageCoordinates::AngularDistance<double>(phaseCentreRA, phaseCentreDec,
				cluster.MeanRA(), cluster.MeanDec()) * (180.0/M_PI);
			
		if(_distanceLimit > 0.0 && distance > _distanceLimit)
		{
			std::cout << "Cluster " << clusterNames[clusterIter] << " is " << distance << " deg away: removing cluster.\n";
			clustersToRemove.insert(clusterIter);
		}
		else if(_clusterFluxLimit > 0.0 && _predicters[clusterIter]->TotalFlux(0) < _clusterFluxLimit)
		{
			std::cout << "Direction for " << clusterNames[clusterIter] << " is below cluster flux limit (" << _predicters[clusterIter]->TotalFlux(0) << " < " << _clusterFluxLimit << "): removing cluster.\n";
			clustersToRemove.insert(clusterIter);
		}
	}
	
	for(std::set<size_t>::const_reverse_iterator i=clustersToRemove.rbegin(); i!=clustersToRemove.rend(); ++i)
	{
		const size_t s = *i;
		_predictionModels.erase(_predictionModels.begin() + s);
		delete _predicters[s];
		_predicters.erase(_predicters.begin() + s);
	}
	std::cout << "Filtered model has " << _predictionModels.size() << " clusters.\n";
	
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
	_channelBlockCount = _bandData.ChannelCount() / _channelBlockSize;
	std::cout << "Will process " << channelCount << " channels in " << _channelBlockCount << " channel groups\n";
	
	_solutionFile.reset(new IonSolutionFile());
	_solutionFile->SetAntennaCount(1);
	_solutionFile->SetChannelBlockCount(_channelBlockCount);
	_solutionFile->SetPolarizationCount(1);
	_solutionFile->SetIntervalCount(_passCount);
	_solutionFile->SetDirectionCount(_predictionModels.size());
	_solutionFile->SetStartFrequency(_bandData.LowestFrequency());
	_solutionFile->SetEndFrequency(_bandData.HighestFrequency());
	_solutionFile->OpenForWriting(solutionFilename);
	
	_model.GetClusterNames(clusterNames);
	for(size_t clusterIter=0; clusterIter!=clusterNames.size(); ++clusterIter)
	{
		SourceGroup cluster;
		_model.GetSourcesInCluster(clusterNames[clusterIter], cluster);
		std::vector<std::string> sourceNames;
		sourceNames.reserve(cluster.SourceCount());
		for(SourceGroup::const_iterator s=cluster.begin(); s!=cluster.end(); ++s)
			sourceNames.push_back(s->Name());
		_solutionFile->WriteClusterMetaInfo(clusterNames[clusterIter], sourceNames);
	}
	
	casa::Array<std::complex<float> > data(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	_dataArrays.resize(_channelBlockCount);
	_weightArrays.resize(_channelBlockCount);
	for(_pass=0; _pass!=_passCount; ++_pass)
	{
		_failedConvergencesPerSource.assign(_predictionModels.size(), 0);
		_failedConvergencesPerChannelGroup.assign(_channelBlockCount, 0);

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
			for(size_t s=0; s!=_predictionModels.size(); ++s)
				_predicters[s]->UpdateBeam(_predictionModels[s]);
		}
		
		std::cout << "Reading (T " << _startTimestep << "-" << _endTimestep << "/" << timestepCount << ", rows " << _curStartRow << '-' << _curEndRow << '/' << ms.nrow() << ")...\n";
		for(size_t channelBlockIndex=0; channelBlockIndex!=_channelBlockCount; ++channelBlockIndex)
		{
			size_t
				channelIndexStart = channelBlockIndex * _bandData.ChannelCount() / _channelBlockCount,
				channelIndexEnd = (channelBlockIndex+1) * _bandData.ChannelCount() / _channelBlockCount,
				channelBlockSize = channelIndexEnd - channelIndexStart;
			_dataArrays[channelBlockIndex] = new VisibilityArray<std::complex<double>, 4>(1, _antennaCount, timestepsInPass * channelBlockSize);
			_weightArrays[channelBlockIndex] = new VisibilityArray<double, 2>(1, _antennaCount, timestepsInPass * channelBlockSize);
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
				
				for(size_t cb=0; cb!=_channelBlockCount; ++cb)
				{
					size_t
						channelIndexStart = cb * _bandData.ChannelCount() / _channelBlockCount,
						channelIndexEnd = (cb+1) * _bandData.ChannelCount() / _channelBlockCount;
					for(size_t channelIndex=channelIndexStart; channelIndex!=channelIndexEnd; ++channelIndex)
					{
						size_t timeIndexOffset = (channelIndex-channelIndexStart)*(_endTimestep-_startTimestep);
						for(size_t p=0; p!=4; ++p)
						{
							if(flagPtr[p]) weightsReadPtr[p] = 0.0;
						}
						
						std::complex<double> *arrPtr = _dataArrays[cb]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep + timeIndexOffset);
						arrPtr[0] = dataPtr[0]; arrPtr[1] = dataPtr[1];
						arrPtr[2] = dataPtr[2]; arrPtr[3] = dataPtr[3];
						double* weightsWritePtr = _weightArrays[cb]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep + timeIndexOffset);
						weightsWritePtr[0] = weightsReadPtr[0];
						weightsWritePtr[1] = weightsReadPtr[3];
						
						weightsReadPtr += 4;
						dataPtr += 4;
						flagPtr += 4;
					}
				}
			}
		}
		
		_progressBar.reset(new ProgressBar("Processing channels"));
		_stats = PeelingStats();
		std::vector<size_t> tasks;
		for(size_t cb=0; cb!=_channelBlockCount; ++cb)
			tasks.push_back(_channelBlockCount - cb - 1);
		std::mutex mutex;
		boost::thread_group threads;
		for(size_t i=0; i!=_cpuCount; ++i)
			threads.add_thread(new boost::thread(&IonPeeler::processingThreadFunction, this, &mutex, &tasks));
		threads.join_all();
		
		// Write back
		for(size_t rowIndex=_curStartRow; rowIndex!=_curEndRow; ++rowIndex)
		{
			const RowData& rowData = _rowData[rowIndex-_curStartRow];
			if(rowData.a1 != rowData.a2)
			{
				timeIndex = rowData.timeIndex;
				std::complex<float> *dataPtr = data.cbegin();
				for(size_t cb=0; cb!=_channelBlockCount; ++cb)
				{
					size_t
						channelIndexStart = cb * _bandData.ChannelCount() / _channelBlockCount,
						channelIndexEnd = (cb+1) * _bandData.ChannelCount() / _channelBlockCount;
					for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
					{
						size_t timeIndexOffset = (ch-channelIndexStart)*(_endTimestep-_startTimestep);
						std::complex<double> *arrPtr = _dataArrays[cb]->ValuePtr(rowData.a1, rowData.a2, timeIndex - _startTimestep + timeIndexOffset);
						dataPtr[0] = arrPtr[0]; dataPtr[1] = arrPtr[1];
						dataPtr[2] = arrPtr[2]; dataPtr[3] = arrPtr[3];
						dataPtr += 4;
					}
				}
				dataColumn.put(rowIndex, data);
			}
		}
		_progressBar.reset();
		
		outputStats(_stats);
		
		for(size_t cb=0; cb!=_channelBlockCount; ++cb)
		{
			if(_failedConvergencesPerChannelGroup[cb] != 0)
			{
				std::cout << "Warning: Solutions for channel group " << cb << " failed to converge " << _failedConvergencesPerChannelGroup[cb] << "x within 100 iterations.\n";
			}
			delete _dataArrays[cb];
			delete _weightArrays[cb];
		}
		
		for(size_t s=0; s!=_predictionModels.size(); ++s)
		{
			if(_failedConvergencesPerSource[s] != 0)
			{
				const ModelSource& firstSource = _predictionModels.front().Source(0);
				std::cout << "Warning: Solutions for " << firstSource.ClusterName() << " failed to converge " << _failedConvergencesPerSource[s] << "x within 100 iterations.\n";
			}
		}
	}
	
	for(size_t s=0; s!=_predictionModels.size(); ++s)
		delete _predicters[s];
}

void IonPeeler::processingThreadFunction(std::mutex* mutex, std::vector<size_t>* tasks)
{
	PeelingStats stats;
	std::unique_lock<std::mutex> lock(*mutex);
	while(!tasks->empty())
	{
		_progressBar->SetProgress(_channelBlockCount - tasks->size(), _channelBlockCount);
		const size_t channel = tasks->back();
		tasks->pop_back();
		lock.unlock();
		
		processChannel(channel, stats);
		
		lock.lock();
	}
	// (mutex must remain locked)
	_stats += stats;
}

void IonPeeler::processChannel(size_t channelIndex, PeelingStats& stats)
{
	positionFitter(channelIndex, stats);
}

// TODO : add channel start and end to fitting Info, and use all channels in fit
int IonPeeler::posMinimizationFunc(const gsl_vector *xvec, void *data, gsl_vector *f)
{
	const FittingInfo &fittingInfo = *reinterpret_cast<FittingInfo*>(data);
	const IonPeeler &ionPeeler = *reinterpret_cast<IonPeeler*>(fittingInfo.ionPeeler);
	const size_t startTimestep = ionPeeler._startTimestep;
	double g = gsl_vector_get(xvec, 0);
	double dl = gsl_vector_get(xvec, 1);
	double dm = gsl_vector_get(xvec, 2);
	
	size_t dataPointIndex = 0;
	const std::complex<double>* modelPtr = &(*fittingInfo.modelData)[0];
	for(size_t ch=fittingInfo.startChannel; ch!=fittingInfo.endChannel; ++ch)
	{
		const double lambda = ionPeeler._bandData.ChannelWavelength(ch);
		size_t timeIndexOffset = (ch-fittingInfo.startChannel)*(ionPeeler._endTimestep-ionPeeler._startTimestep);
		for(size_t row=ionPeeler._curStartRow; row!=ionPeeler._curEndRow; ++row)
		{
			const RowData& rowData = ionPeeler._rowData[row - ionPeeler._curStartRow];
			if(rowData.a1 != rowData.a2)
			{
				const std::complex<double>* dataPtr = ionPeeler._dataArrays[fittingInfo.channelBlockIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimestep + timeIndexOffset);
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
	}
		
	return GSL_SUCCESS;
}
int IonPeeler::posMinimizationFuncDeriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
{
	const FittingInfo &fittingInfo = *reinterpret_cast<FittingInfo*>(data);
	const IonPeeler &ionPeeler = *reinterpret_cast<IonPeeler*>(fittingInfo.ionPeeler);
	const size_t startTimestep = ionPeeler._startTimestep;
	
	double g = gsl_vector_get(xvec, 0);
	double dl = gsl_vector_get(xvec, 1);
	double dm = gsl_vector_get(xvec, 2);
	
	size_t dataPointIndex = 0;
	const std::complex<double>* modelPtr = &(*fittingInfo.modelData)[0];
	for(size_t ch=fittingInfo.startChannel; ch!=fittingInfo.endChannel; ++ch)
	{
		const double lambda = ionPeeler._bandData.ChannelWavelength(ch);
		size_t timeIndexOffset = (ch-fittingInfo.startChannel)*(ionPeeler._endTimestep-ionPeeler._startTimestep);
		for(size_t row=ionPeeler._curStartRow; row!=ionPeeler._curEndRow; ++row)
		{
			const RowData& rowData = ionPeeler._rowData[row - ionPeeler._curStartRow];
			if(rowData.a1 != rowData.a2)
			{
				std::complex<double> dfdg = 0.0, dfddl = 0.0, dfddm = 0.0;
				const std::complex<double>* dataPtr = ionPeeler._dataArrays[fittingInfo.channelBlockIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - startTimestep + timeIndexOffset);
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
	}
	return GSL_SUCCESS;
}
int IonPeeler::posMinimizationFuncBoth(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	posMinimizationFunc(x, data, f);
	posMinimizationFuncDeriv(x, data, J);
	return GSL_SUCCESS;
}

void IonPeeler::positionFitter(size_t channelBlockIndex, PeelingStats& stats)
{
	size_t
		channelIndexStart = channelBlockIndex * _bandData.ChannelCount() / _channelBlockCount,
		channelIndexEnd = (channelBlockIndex+1) * _bandData.ChannelCount() / _channelBlockCount,
		curChannelBlockSize = channelIndexEnd - channelIndexStart;
	
	FittingInfo fInfo;
	fInfo.ionPeeler = this;
	fInfo.channelBlockIndex = channelBlockIndex;
	fInfo.startChannel = channelIndexStart;
	fInfo.endChannel = channelIndexEnd;
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	size_t dataValueCount = 0;
	for(size_t row=_curStartRow; row!=_curEndRow; ++row) 
	{
		const RowData& rowData = _rowData[row - _curStartRow];
		if(rowData.a1 != rowData.a2) ++dataValueCount;
	}
	dataValueCount *= 2 /* For real and imaginary */ * (channelIndexEnd-channelIndexStart);
	gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc(T, dataValueCount, 3);
	
	ao::uvector<double>
		_solutionsG(_predictionModels.size()),
		_solutionsDL(_predictionModels.size()),
		_solutionsDM(_predictionModels.size());
	
	// The array containing for each source, the model data per cross-correlated row: modelData[s][row*4+p]
	std::vector<ao::uvector<std::complex<double>>> modelDataPerSource(_predictionModels.size());
	
	for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
	{
		ao::uvector<std::complex<double>>& modelData = modelDataPerSource[sourceIndex];
		modelData.resize((_curEndRow-_curStartRow)*4*curChannelBlockSize);
	}
		
	fInfo.lambda = 0.0;
	for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
	{
		const double lambda = _bandData.ChannelWavelength(ch);
		fInfo.lambda += lambda;
	}
	fInfo.lambda /= curChannelBlockSize;
		
	// Predict model for this source
	for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
	{
		ao::uvector<std::complex<double>>& modelData = modelDataPerSource[sourceIndex];
		ao::uvector<std::complex<double>>::iterator modelPtr = modelData.begin();
		for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
		{
			const double lambda = _bandData.ChannelWavelength(ch);
			for(size_t row=_curStartRow; row!=_curEndRow; ++row)
			{
				const RowData& rowData = _rowData[row - _curStartRow];
				if(rowData.a1 != rowData.a2)
				{
					const double uOverL = rowData.u/lambda, vOverL = rowData.v/lambda, wOverL = rowData.w/lambda;
					_predicters[sourceIndex]->Predict4(&*modelPtr, _predictionModels[sourceIndex], uOverL, vOverL, wOverL, ch, rowData.a1, rowData.a2);
					modelPtr += 4;
				}
			}
		}
	}
	fInfo.lambda /= curChannelBlockSize;
	
	for(size_t fitIteration=0; fitIteration!=_fitIterationCount; ++fitIteration)
	{
		for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
		{
			ao::uvector<std::complex<double>>& modelData = modelDataPerSource[sourceIndex];
			fInfo.modelData = &modelData;
			double
				&g = _solutionsG[sourceIndex],
				&dl = _solutionsDL[sourceIndex],
				&dm = _solutionsDM[sourceIndex];
			
			/*if(channelIndex == 5)
			{
				std::cout << "Start of fit (" << sourceIndex << "): " <<
					g << ',' << dl << ',' << dm << '\n';
			}*/
				
			// Add back if it was subtracted
			if(fitIteration != 0)
			{
				ao::uvector<std::complex<double>>::iterator modelPtr = modelData.begin();
				for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
				{
					const double lambda = _bandData.ChannelWavelength(ch);
					size_t timeIndexOffset = (ch-channelIndexStart)*(_endTimestep-_startTimestep);
					for(size_t row=_curStartRow; row!=_curEndRow; ++row)
					{
						const RowData& rowData = _rowData[row - _curStartRow];
						if(rowData.a1 != rowData.a2)
						{
							std::complex<double>* dataPtr = _dataArrays[channelBlockIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - _startTimestep + timeIndexOffset);
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
								dataPtr[p] += corModel;
							}
							modelPtr += 4;
						}
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
			
			if(iter == 100)
			{
				_failedConvergencesPerSource[sourceIndex]++;
				_failedConvergencesPerChannelGroup[channelBlockIndex]++;
			}
			
			g = gsl_vector_get (solver->x, 0);
			dl = gsl_vector_get (solver->x, 1);
			dm = gsl_vector_get (solver->x, 2);
			if(sourceIndex <= 4)
				std::cout << "Solution " << fitIteration << ", channels " << channelIndexStart << '-' << channelIndexEnd << " for " << _predictionModels[sourceIndex].Source(0).ClusterName() << ": " << g << ',' << dl << ',' << dm << '\n';
			//std::cout << gsl_strerror(status) << "\n";
			//std::cout << (status==GSL_CONTINUE ? "CONTINUE " : "X ") << iter << ", g=" << g << ",dl=" << dl << ",dm=-" << dm << '\n';
			/*if(channelIndex == 5)
			{
				std::cout << "End of fit: (" << sourceIndex << ")" <<
					g << ',' << dl << ',' << dm << '\n';
			}*/
				
			
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
			ao::uvector<std::complex<double>>::iterator modelPtr = modelData.begin();
			for(size_t ch=channelIndexStart; ch!=channelIndexEnd; ++ch)
			{
				const double lambda = _bandData.ChannelWavelength(ch);
				size_t timeIndexOffset = (ch-channelIndexStart)*(_endTimestep-_startTimestep);
				for(size_t row=_curStartRow; row!=_curEndRow; ++row)
				{
					const RowData& rowData = _rowData[row - _curStartRow];
					if(rowData.a1 != rowData.a2)
					{
						std::complex<double>* dataPtr = _dataArrays[channelBlockIndex]->ValuePtr(rowData.a1, rowData.a2, rowData.timeIndex - _startTimestep + timeIndexOffset);
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
	}
	gsl_multifit_fdfsolver_free(solver);
	
	std::vector<IonSolutionFile::Solution> solutions(_predictionModels.size());
	for(size_t sourceIndex=0; sourceIndex!=_predictionModels.size(); ++sourceIndex)
	{
		IonSolutionFile::Solution &solution = solutions[sourceIndex];
		solution.gain = _solutionsG[sourceIndex];
		solution.dl = _solutionsDL[sourceIndex];
		solution.dm = _solutionsDM[sourceIndex];
	}	
	_solutionFile->WriteChannelBlock(solutions.data(), _pass, channelBlockIndex, 0);
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
