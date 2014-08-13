#include "calibrationmethod.h"
#include "matrix2x2.h"

#include <iostream>
#include <cmath>
#include <sstream>
#include <limits>

CalibrationMethod::CalibrationMethod(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
	_data(nChannels, nAntenna, nTimesteps),
	_model(nChannels, nAntenna, nTimesteps),
	_weights(nChannels, nAntenna, nTimesteps),
	_weightSums(1, nAntenna, 1),
	_jonesSolutions(nAntenna * 4 * nChannels, 0.0),
	_nChannels(nChannels),
	_nAntenna(nAntenna),
	_nTimesteps(nTimesteps),
	_onlySolveDiag(false),
	_onlySolveScalar(false),
	_onlySolveRotation(false)
{
	InitSolutionsToUnity();
}

void CalibrationMethod::InitSolutionsToUnity()
{
	std::complex<double> *_jonesPtr = &_jonesSolutions[0];
	for(size_t i=0; i!=_nAntenna * _nChannels; ++i)
	{
		*_jonesPtr = std::complex<double>(1.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(0.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(0.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(1.0, 0.0); ++_jonesPtr;
	}
}

void CalibrationMethod::InitSolutionsToNaN()
{
	std::complex<double> *_jonesPtr = &_jonesSolutions[0];
	const std::complex<double> nan(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	for(size_t i=0; i!=_nAntenna * _nChannels; ++i)
	{
		*_jonesPtr = nan; ++_jonesPtr;
		*_jonesPtr = nan; ++_jonesPtr;
		*_jonesPtr = nan; ++_jonesPtr;
		*_jonesPtr = nan; ++_jonesPtr;
	}
}

void CalibrationMethod::InitSolutions(const CalibrationMethod &source)
{
	_jonesSolutions = source._jonesSolutions;
	ao::uvector<std::complex<double>>::iterator ptr = _jonesSolutions.begin();
	for(size_t i=0; i!=_nChannels*_nAntenna; ++i)
	{
		bool isFlagged = false;
		for(size_t p=0; p!=4; ++p)
		{
			if(std::isfinite(ptr->real()) || std::isfinite(ptr->imag()))
				isFlagged = true;
			++ptr;
		}
		if(isFlagged)
		{
			*(ptr-4) = 1.0; *(ptr-3) = 0.0;
			*(ptr-2) = 0.0; *(ptr-1) = 1.0;
		}
	}
}

template<typename DataFloatType>
void CalibrationMethod::AddData(const std::complex<DataFloatType>* data, const float* weights, const std::complex<double>* predictedValues, size_t antenna1, size_t antenna2, size_t timestep)
{
	std::complex<double> *destDataPtr = _data.ValuePtr(antenna1, antenna2, timestep);
	std::complex<double> *destModelPtr = _model.ValuePtr(antenna1, antenna2, timestep);
	double *destWeightPtr = _weights.ValuePtr(antenna1, antenna2, timestep);
	
	const std::complex<DataFloatType>* dataEndPtr = data + _nChannels * 4;
	
	while(data != dataEndPtr)
	{
		double minWeight = *weights;
		for(size_t i=0; i!=4; ++i)
		{
			if(std::isfinite(data->real()) && std::isfinite(data->imag()))
			{
				*destDataPtr = std::complex<double>(*data);
				*destModelPtr = *predictedValues;
				if(*weights < minWeight)
					minWeight = *weights;
			}
			else {
				*destDataPtr = std::complex<double>(0.0, 0.0);
				*destModelPtr = std::complex<double>(0.0, 0.0);
				minWeight = 0.0;
			}
			
			++data;
			++weights;
			++predictedValues;
			++destDataPtr;
			++destModelPtr;
		}
		
		*destWeightPtr = minWeight;
		++destWeightPtr;
	}
}

template
void CalibrationMethod::AddData(const std::complex<float> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2, size_t timestep);

template
void CalibrationMethod::AddData(const std::complex<double> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2, size_t timestep);


void CalibrationMethod::applyWeightsToData()
{
	_weightSums.SetAll(0.0);
	
	for(size_t timestep=0; timestep!=_nTimesteps; ++timestep)
	{
		for(size_t antenna2 = 0; antenna2!=_nAntenna; ++antenna2)
		{
			for(size_t antenna1 = 0; antenna1!=antenna2; ++antenna1)
			{
				std::complex<double> *dataPtr = _data.ValuePtr(antenna1, antenna2, timestep);
				std::complex<double> *modelPtr = _model.ValuePtr(antenna1, antenna2, timestep);
				double *weightPtr = _weights.ValuePtr(antenna1, antenna2, timestep);
				double &weightSum = *_weightSums.ValuePtr(antenna1, antenna2, 0);
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					const double w = *weightPtr;
					weightSum += w;
					for(size_t p=0; p!=4; ++p)
					{
						*dataPtr *= w;
						*modelPtr *= w;
						++dataPtr;
						++modelPtr;
					}
					++weightPtr;
				}
			}
		}
	}
}

double CalibrationMethod::totalDistance(size_t antenna)
{
	double distance = 0.0, weightSum = 0.0;
	std::complex<double> xxFlux = 0.0, modelXX = 0.0;
	size_t xxFluxCount = 0;
	for(size_t timestep=0; timestep!=_nTimesteps; ++timestep)
	{
		for(size_t antenna1 = 0; antenna1!=_nAntenna; ++antenna1)
		{
			if(antenna1 != antenna)
			{
				const bool isConjTranspose = antenna > antenna1;
				const std::complex<double> *data = _data.ValuePtr(antenna, antenna1, timestep);
				const std::complex<double> *model = _model.ValuePtr(antenna, antenna1, timestep);
				const std::complex<double> *jones1 = &_jonesSolutions[antenna * _nChannels * 4];
				const std::complex<double> *jones2 = &_jonesSolutions[antenna1 * _nChannels * 4];
				
				std::complex<double> j1TimesD[4], j1DJ2[4], distances[4];
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					if(isConjTranspose)
					{
						j1TimesD[0] = jones1[0] * std::conj(data[0]) + jones1[1] * std::conj(data[1]);
						j1TimesD[1] = jones1[0] * std::conj(data[2]) + jones1[1] * std::conj(data[3]);
						j1TimesD[2] = jones1[2] * std::conj(data[0]) + jones1[3] * std::conj(data[1]);
						j1TimesD[3] = jones1[2] * std::conj(data[2]) + jones1[3] * std::conj(data[3]);
						
						j1DJ2[0] = j1TimesD[0] * std::conj(jones2[0/* (0^H) */]) + j1TimesD[1] * std::conj(jones2[1/* (2^H) */]);
						j1DJ2[1] = j1TimesD[0] * std::conj(jones2[2/* (1^H) */]) + j1TimesD[1] * std::conj(jones2[3/* (3^H) */]);
						j1DJ2[2] = j1TimesD[2] * std::conj(jones2[0/* (0^H) */]) + j1TimesD[3] * std::conj(jones2[1/* (2^H) */]);
						j1DJ2[3] = j1TimesD[2] * std::conj(jones2[2/* (1^H) */]) + j1TimesD[3] * std::conj(jones2[3/* (3^H) */]);
						
						xxFlux += j1DJ2[0]; xxFluxCount++; modelXX += model[0];
						distances[0] = std::conj(model[0]) - j1DJ2[0];
						distances[1] = std::conj(model[2]) - j1DJ2[1],
						distances[2] = std::conj(model[1]) - j1DJ2[2];
						distances[3] = std::conj(model[3]) - j1DJ2[3];
					} else {
						j1TimesD[0] = jones1[0] * data[0] + jones1[1] * data[2];
						j1TimesD[1] = jones1[0] * data[1] + jones1[1] * data[3];
						j1TimesD[2] = jones1[2] * data[0] + jones1[3] * data[2];
						j1TimesD[3] = jones1[2] * data[1] + jones1[3] * data[3];
						
						j1DJ2[0] = j1TimesD[0] * std::conj(jones2[0/* (0^H) */]) + j1TimesD[1] * std::conj(jones2[1/* (2^H) */]);
						j1DJ2[1] = j1TimesD[0] * std::conj(jones2[2/* (1^H) */]) + j1TimesD[1] * std::conj(jones2[3/* (3^H) */]);
						j1DJ2[2] = j1TimesD[2] * std::conj(jones2[0/* (0^H) */]) + j1TimesD[3] * std::conj(jones2[1/* (2^H) */]);
						j1DJ2[3] = j1TimesD[2] * std::conj(jones2[2/* (1^H) */]) + j1TimesD[3] * std::conj(jones2[3/* (3^H) */]);
						
						xxFlux += j1DJ2[0]; xxFluxCount++; modelXX += model[0];
						distances[0] = model[0] - j1DJ2[0];
						distances[1] = model[1] - j1DJ2[1],
						distances[2] = model[2] - j1DJ2[2];
						distances[3] = model[3] - j1DJ2[3];
					}
					
					distance +=
						std::norm(distances[0]) + std::norm(distances[1]) +
						std::norm(distances[2]) + std::norm(distances[3]);
					
					data += 4;
					model += 4;
					jones1 += 4;
					jones2 += 4;
				}
				
				weightSum += *_weightSums.ValuePtr(antenna, antenna1, 0);
			}
		}
	}
	std::cout << "XX flux for ant " << antenna << ": " << (xxFlux / double(xxFluxCount)) << " model=" << (modelXX / double(xxFluxCount)) << '\n';
	if(weightSum == 0.0)
		return 0.0;
	else
		return sqrt(distance / (weightSum * 4));
}

double CalibrationMethod::totalDistance()
{
	double distance = 0.0;
	for(size_t antenna = 0; antenna!=_nAntenna; ++antenna)
	{
		distance += totalDistance(antenna);
	}
	return distance;
}

void CalibrationMethod::reportDistances()
{
	double sumDistance = 0.0;
	std::cout << "Distances: ";
	for(size_t ant=0; ant!=_nAntenna; ++ant)
	{
		double thisDistance = totalDistance(ant);
		std::cout << "a" << ant << ": " << thisDistance << '\t' << std::flush;
		sumDistance += thisDistance;
	}
	std::cout << "\nTotal average: " << (sumDistance / _nAntenna) << '\n';
}

void CalibrationMethod::Execute(double& precisionLimit, size_t& nIter)
{
	bool continueIterating;
	size_t iterationNumber = 0;
	
	_weightSums.SetAll(_nChannels * _nTimesteps);
	//reportDistances();
	
	//std::cout << "Weighting data.\n";
	applyWeightsToData();
	
	double globalChangeSizes[4] = {0.0,0.0,0.0,0.0};
	double stepsize = 0.25;
	ao::uvector<bool> antennaResults(_nAntenna * _nChannels, false);
	do
	{
		++iterationNumber;
		//std::cout << "Iteration " << iterationNumber << '\n';
		
		//reportDistances();
		
		ao::uvector<std::complex<double> > nextJones(_jonesSolutions);
		
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			calculateNextIter(ant, &nextJones[ant*4*_nChannels], &antennaResults[ant*_nChannels]);
			if(_onlySolveScalar) {
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					nextJones[(ant*_nChannels+ch)*4 + 0] = (nextJones[(ant*_nChannels+ch)*4 + 0] + nextJones[(ant*_nChannels+ch)*4 + 3]) * 0.5;
					nextJones[(ant*_nChannels+ch)*4 + 1] = 0;
					nextJones[(ant*_nChannels+ch)*4 + 2] = 0;
					nextJones[(ant*_nChannels+ch)*4 + 3] = nextJones[(ant*_nChannels+ch)*4 + 0];
				}
			}
			else if(_onlySolveDiag) {
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					nextJones[(ant*_nChannels+ch)*4 + 1] = 0;
					nextJones[(ant*_nChannels+ch)*4 + 2] = 0;
				}
			}
		}
		
		std::complex<double> *jonesPtr = &_jonesSolutions[0];
		std::complex<double> *nextJonesPtr = &nextJones[0];
		std::vector<double> changeSizes(_nAntenna*4);
		// TODO stepsize based on something
		//stepsize *= 0.99;
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			for(size_t ch=0; ch!=_nChannels; ++ch)
			{
				if(_onlySolveRotation) {
					double alpha = Matrix2x2::RotationAngle(nextJonesPtr);
					std::complex<double> temp[4];
					Matrix2x2::RotationMatrix(temp, alpha);
					for(size_t p=0; p!=4; ++p)
					{
						*jonesPtr = *jonesPtr * (1.0-stepsize) + temp[p] * stepsize;
						changeSizes[p+ant*4] += std::norm(*jonesPtr - temp[p]);
						++jonesPtr;
					}
					nextJonesPtr += 4;
				}
				else {
					for(size_t p=0; p!=4; ++p)
					{
						changeSizes[p+ant*4] += std::norm(*jonesPtr - *nextJonesPtr);
						
						*jonesPtr = *jonesPtr * (1.0-stepsize) + *nextJonesPtr * stepsize;
						
						++jonesPtr;
						++nextJonesPtr;
					}
				}
				

				/*if(ant==1 && (ch==15 || ch==14 || ch==13))
				{
					std::cout << "Current value of Jones matrix for ant 1, ch " << ch << ":\n"
					<< matrixToString(jonesPtr-4);
				}*/
			}
		}
		
		continueIterating = false;
		for(size_t p=0; p!=4; ++p)
			globalChangeSizes[p] = 0.0;
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			for(size_t p=0; p!=4; ++p)
			{
				changeSizes[p+ant*4] /= _nChannels;
				globalChangeSizes[p] += changeSizes[p+ant*4];
			}
		}
		for(size_t p=0; p!=4; ++p)
		{
			globalChangeSizes[p] /= _nAntenna;
			if(globalChangeSizes[p] > precisionLimit)
				continueIterating = true;
		}
		
		//std::cout << "Average change to Jones solutions: \n"
		//	" (" << globalChangeSizes[0] << " " << globalChangeSizes[1] << ")\n"
		//	" (" << globalChangeSizes[2] << " " << globalChangeSizes[3] << ")\n";
		
	} while(continueIterating && iterationNumber<nIter);

	//reportDistances();
	nIter = iterationNumber;
	precisionLimit = std::max(
		std::max(globalChangeSizes[0], globalChangeSizes[1]),
		std::max(globalChangeSizes[2], globalChangeSizes[3]));
	
	std::complex<double> *jonesPtr = &_jonesSolutions[0];
	const std::complex<double> nan(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
	for(size_t ant=0; ant!=_nAntenna; ++ant)
	{
		for(size_t ch=0; ch!=_nChannels; ++ch)
		{
			if(!antennaResults[ant*_nChannels + ch])
			{
				for(size_t p=0; p!=4; ++p)
				{
					*jonesPtr = nan;
					++jonesPtr;
				}
			}
			else {
				jonesPtr += 4;
			}
		}
	}
}

std::string CalibrationMethod::MatrixToString(const std::complex<double> *matrix)
{
	std::stringstream s;
	double s1, s2;
	Matrix2x2::SingularValues(matrix, s1, s2);
	s <<
		" (" << matrix[0] << " " << matrix[1] << ")\tamplitudes=(" << abs(matrix[0]) << ' ' << abs(matrix[1]) << ") SV=\t" << s1 << "\n"
		" (" << matrix[2] << " " << matrix[3] << ")\t           (" << abs(matrix[2]) << ' ' << abs(matrix[3]) << ")    \t" << s2 << '\n';
	return s.str();
}

void CalibrationMethod::calculateNextIter(size_t ant, std::complex<double> *nextJones, bool* antennaResults)
{
	// Calculate the next Jones estimate for the given ant, based on
	// the previous Jones estimates.
	
	// Calculate:
	// JONES[ant] = SUM_{k!=ant in Nant} MODEL[ant,k] JONES[k] DATA^H[ant,k] times
	//            ( SUM_{k!=ant in Nant} DATA[ant,k] JONES^H[k] JONES[k] DATA^H[ant,k] )^-1
	// (From Mitchel et al., 2008)
	
	std::vector<std::complex<double> > solutions(4 * _nChannels);
	
	std::vector<std::complex<double> > rTerm(4 * _nChannels);
	for(size_t t=0; t!=_nTimesteps; ++t)
	{
		for(size_t k=0; k!=_nAntenna; ++k)
		{
			if(k != ant)
			{
				const bool isConjTranspose = ant > k;
				std::complex<double> *dataPtr = _data.ValuePtr(ant, k, t);
				std::complex<double> *modelPtr = _model.ValuePtr(ant, k, t);
				std::complex<double> *jonesPtr = &_jonesSolutions[k * 4 * _nChannels];
				std::complex<double> *rTermPtr = &rTerm[0];
				std::complex<double> *nextJonesPtr = &solutions[0];
			
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					if(isConjTranspose)
					{
						// sum(D^H J M) sum(M^H J^H J M)
						std::complex<double> jTimesHermM[4];
						
						Matrix2x2::ATimesB(jTimesHermM, jonesPtr, modelPtr);
						
						Matrix2x2::PlusHermATimesB(nextJonesPtr, dataPtr, jTimesHermM);
						
						std::complex<double> mTimesHermJ[4];
						Matrix2x2::HermATimesHermB(mTimesHermJ, modelPtr, jonesPtr);
						
						Matrix2x2::PlusATimesB(rTermPtr, mTimesHermJ, jTimesHermM);
						
					} else { // non-herm conjugate case
						// sum(D J M^H) sum(M J^H J M^H)
						std::complex<double> jTimesHermM[4];
						
						Matrix2x2::ATimesHermB(jTimesHermM, jonesPtr, modelPtr);
						
						Matrix2x2::PlusATimesB(nextJonesPtr, dataPtr, jTimesHermM);
						
						std::complex<double> mTimesHermJ[4];
						Matrix2x2::ATimesHermB(mTimesHermJ, modelPtr, jonesPtr);
						
						Matrix2x2::PlusATimesB(rTermPtr, mTimesHermJ, jTimesHermM);
					}
					
					// Move to next channel
					nextJonesPtr += 4;
					dataPtr += 4;
					modelPtr += 4;
					jonesPtr += 4;
					rTermPtr += 4;
				}
			}
		}
	}
	
	for(size_t ch=0; ch!=_nChannels; ++ch)
	{
		if(Matrix2x2::MultiplyWithInverse(&solutions[ch * 4], &rTerm[ch * 4]))
		{
			antennaResults[ch] = true;
			for(size_t i=0; i!=4; ++i)
				nextJones[ch * 4 + i] = solutions[ch * 4 + i];
		}
		else antennaResults[ch] = false;
	}
}

void CalibrationMethod::SolutionSingularValue(size_t antenna, size_t channel, double &s1, double &s2) const
{
	Matrix2x2::SingularValues(&_jonesSolutions[(antenna * _nChannels + channel) * 4], s1, s2);
}
