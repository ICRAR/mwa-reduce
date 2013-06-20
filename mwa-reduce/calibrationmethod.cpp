#include "calibrationmethod.h"

#include <iostream>
#include <cmath>
#include <sstream>

CalibrationMethod::CalibrationMethod(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
	_data(nChannels, nAntenna, nTimesteps),
	_model(nChannels, nAntenna, nTimesteps),
	_weights(nChannels, nAntenna, nTimesteps),
	_weightSums(1, nAntenna, 1),
	_jonesSolutions(nAntenna * 4 * nChannels),
	_nChannels(nChannels),
	_nAntenna(nAntenna),
	_nTimesteps(nTimesteps)
{
	std::complex<double> *_jonesPtr = &_jonesSolutions[0];
	for(size_t i=0; i!=nAntenna * nChannels; ++i)
	{
		*_jonesPtr = std::complex<double>(1.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(0.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(0.0, 0.0); ++_jonesPtr;
		*_jonesPtr = std::complex<double>(1.0, 0.0); ++_jonesPtr;
	}
}

void CalibrationMethod::AddData(const std::complex<float>* data, const float* weights, const std::complex<double>* predictedValues, size_t antenna1, size_t antenna2, size_t timestep)
{
	std::complex<double> *destDataPtr = _data.ValuePtr(antenna1, antenna2, timestep);
	std::complex<double> *destModelPtr = _model.ValuePtr(antenna1, antenna2, timestep);
	double *destWeightPtr = _weights.ValuePtr(antenna1, antenna2, timestep);
	
	const std::complex<float>* dataEndPtr = data + _nChannels * 4;
	
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
				*destModelPtr = *predictedValues;
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

void CalibrationMethod::Execute(double precisionLimit, size_t nIter)
{
	bool continueIterating;
	size_t iterationNumber = 0;
	
	_weightSums.SetAll(_nChannels * _nTimesteps);
	//reportDistances();
	
	std::cout << "Weighting data.\n";
	applyWeightsToData();
	
	do
	{
		++iterationNumber;
		std::cout << "Iteration " << iterationNumber << '\n';
		
		//reportDistances();
		
		std::vector<std::complex<double> > nextJones(_jonesSolutions);
		
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			calculateNextIter(ant, &nextJones[ant*4*_nChannels]);
		}
		
		std::complex<double> *jonesPtr = &_jonesSolutions[0];
		std::complex<double> *nextJonesPtr = &nextJones[0];
		std::vector<double> changeSizes(_nAntenna*4);
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			for(size_t ch=0; ch!=_nChannels; ++ch)
			{
				for(size_t p=0; p!=4; ++p)
				{
					changeSizes[p+ant*4] += std::norm(*jonesPtr - *nextJonesPtr);
					
					// TODO stepsize based on something
					const double STEPSIZE = 0.25;
					*jonesPtr = *jonesPtr * (1.0-STEPSIZE) + *nextJonesPtr * STEPSIZE;
					
					++jonesPtr;
					++nextJonesPtr;
				}
				
				if(ant==1 && (ch==15 || ch==14 || ch==13))
				{
					std::cout << "Current value of Jones matrix for ant 1, ch " << ch << ":\n"
					<< matrixToString(jonesPtr-4);
				}
			}
		}
		
		continueIterating = false;
		double globalChangeSizes[4] = { 0.0, 0.0, 0.0, 0.0 };
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
		
		std::cout << "Average change to Jones solutions: \n"
			" (" << globalChangeSizes[0] << " " << globalChangeSizes[1] << ")\n"
			" (" << globalChangeSizes[2] << " " << globalChangeSizes[3] << ")\n";
		
	} while(continueIterating && iterationNumber<nIter);
}

std::string CalibrationMethod::matrixToString(std::complex<double> *matrix)
{
	std::stringstream s;
	double s1, s2;
	singularValues2x2(matrix, s1, s2);
	s <<
		" (" << matrix[0] << " " << matrix[1] << ")\tamplitudes=(" << abs(matrix[0]) << ' ' << abs(matrix[1]) << ") SV=\t" << s1 << "\n"
		" (" << matrix[2] << " " << matrix[3] << ")\t           (" << abs(matrix[2]) << ' ' << abs(matrix[3]) << ")    \t" << s2 << '\n';
	return s.str();
}

void CalibrationMethod::calculateNextIter(size_t ant, std::complex<double> *nextJones)
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
						std::complex<double> jTimesHermD[4] = {
							(jonesPtr[0] * dataPtr[0] + jonesPtr[1] * dataPtr[2]),
							(jonesPtr[0] * dataPtr[1] + jonesPtr[1] * dataPtr[3]),
							(jonesPtr[2] * dataPtr[0] + jonesPtr[3] * dataPtr[2]),
							(jonesPtr[2] * dataPtr[1] + jonesPtr[3] * dataPtr[3])
						};
						
						nextJonesPtr[0] += std::conj(modelPtr[0]) * jTimesHermD[0] + std::conj(modelPtr[2]) * jTimesHermD[2];
						nextJonesPtr[1] += std::conj(modelPtr[0]) * jTimesHermD[1] + std::conj(modelPtr[2]) * jTimesHermD[3];
						nextJonesPtr[2] += std::conj(modelPtr[1]) * jTimesHermD[0] + std::conj(modelPtr[3]) * jTimesHermD[2];
						nextJonesPtr[3] += std::conj(modelPtr[1]) * jTimesHermD[1] + std::conj(modelPtr[3]) * jTimesHermD[3];
						
						std::complex<double> dTimesHermJ[4] = {
							(std::conj(dataPtr[0]) * std::conj(jonesPtr[0/* (0^H) */]) + std::conj(dataPtr[2]) * std::conj(jonesPtr[1/* (2^H) */])),
							(std::conj(dataPtr[0]) * std::conj(jonesPtr[2/* (1^H) */]) + std::conj(dataPtr[2]) * std::conj(jonesPtr[3/* (3^H) */])),
							(std::conj(dataPtr[1]) * std::conj(jonesPtr[0/* (0^H) */]) + std::conj(dataPtr[3]) * std::conj(jonesPtr[1/* (2^H) */])),
							(std::conj(dataPtr[1]) * std::conj(jonesPtr[2/* (1^H) */]) + std::conj(dataPtr[3]) * std::conj(jonesPtr[3/* (3^H) */]))
						};
						
						rTermPtr[0] += dTimesHermJ[0] * jTimesHermD[0] + dTimesHermJ[1] * jTimesHermD[2];
						rTermPtr[1] += dTimesHermJ[0] * jTimesHermD[1] + dTimesHermJ[1] * jTimesHermD[3];
						rTermPtr[2] += dTimesHermJ[2] * jTimesHermD[0] + dTimesHermJ[3] * jTimesHermD[2];
						rTermPtr[3] += dTimesHermJ[2] * jTimesHermD[1] + dTimesHermJ[3] * jTimesHermD[3];
						
					} else { // non-herm conjugate case
						std::complex<double> jTimesHermD[4] = {
							(jonesPtr[0] * std::conj(dataPtr[0/* (0^H) */]) + jonesPtr[1] * std::conj(dataPtr[1/* (2^H) */])),
							(jonesPtr[0] * std::conj(dataPtr[2/* (1^H) */]) + jonesPtr[1] * std::conj(dataPtr[3/* (3^H) */])),
							(jonesPtr[2] * std::conj(dataPtr[0/* (0^H) */]) + jonesPtr[3] * std::conj(dataPtr[1/* (2^H) */])),
							(jonesPtr[2] * std::conj(dataPtr[2/* (1^H) */]) + jonesPtr[3] * std::conj(dataPtr[3/* (3^H) */]))
						};
						
						nextJonesPtr[0] += modelPtr[0] * jTimesHermD[0] + modelPtr[1] * jTimesHermD[2];
						nextJonesPtr[1] += modelPtr[0] * jTimesHermD[1] + modelPtr[1] * jTimesHermD[3];
						nextJonesPtr[2] += modelPtr[2] * jTimesHermD[0] + modelPtr[3] * jTimesHermD[2];
						nextJonesPtr[3] += modelPtr[2] * jTimesHermD[1] + modelPtr[3] * jTimesHermD[3];
						
						std::complex<double> dTimesHermJ[4] = {
							(dataPtr[0] * std::conj(jonesPtr[0/* (0^H) */]) + dataPtr[1] * std::conj(jonesPtr[1/* (2^H) */])),
							(dataPtr[0] * std::conj(jonesPtr[2/* (1^H) */]) + dataPtr[1] * std::conj(jonesPtr[3/* (3^H) */])),
							(dataPtr[2] * std::conj(jonesPtr[0/* (0^H) */]) + dataPtr[3] * std::conj(jonesPtr[1/* (2^H) */])),
							(dataPtr[2] * std::conj(jonesPtr[2/* (1^H) */]) + dataPtr[3] * std::conj(jonesPtr[3/* (3^H) */]))
						};
						
						rTermPtr[0] += dTimesHermJ[0] * jTimesHermD[0] + dTimesHermJ[1] * jTimesHermD[2];
						rTermPtr[1] += dTimesHermJ[0] * jTimesHermD[1] + dTimesHermJ[1] * jTimesHermD[3];
						rTermPtr[2] += dTimesHermJ[2] * jTimesHermD[0] + dTimesHermJ[3] * jTimesHermD[2];
						rTermPtr[3] += dTimesHermJ[2] * jTimesHermD[1] + dTimesHermJ[3] * jTimesHermD[3];
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
		if(multiplyWithInverse2x2(&solutions[ch * 4], &rTerm[ch * 4]))
		{
			for(size_t i=0; i!=4; ++i)
				nextJones[ch * 4 + i] = solutions[ch * 4 + i];
		}
	}
}

bool CalibrationMethod::multiplyWithInverse2x2(std::complex<double>* lhs, const std::complex<double>* rhs)
{
	std::complex<double> d = ((rhs[0]*rhs[3]) - (rhs[1]*rhs[2]));
	if(d == 0.0) return false;
	std::complex<double> oneOverDeterminant = 1.0 / d;
	std::complex<double> temp[4];
	temp[0] = rhs[3] * oneOverDeterminant;
	temp[1] = -rhs[1] * oneOverDeterminant;
	temp[2] = -rhs[2] * oneOverDeterminant;
	temp[3] = rhs[0] * oneOverDeterminant;
	
	std::complex<double> temp2 = lhs[0];
	lhs[0] = lhs[0] * temp[0] + lhs[1] * temp[2];
	lhs[1] =  temp2 * temp[1] + lhs[1] * temp[3];
	
	temp2 = lhs[2];
	lhs[2] = lhs[2] * temp[0] + lhs[3] * temp[2];
	lhs[3] = temp2 * temp[1] + lhs[3] * temp[3];
	return true;
}

void CalibrationMethod::singularValues2x2(const std::complex<double>* matrix, double &e1, double &e2)
{
	// This is not the ultimate fastest method, since we
	// don't need to calculate the imaginary values of b,c at all.
  // Calculate M M^H
	std::complex<double> temp[4] = {
		matrix[0] * std::conj(matrix[0]) + matrix[1] * std::conj(matrix[1]),
		matrix[0] * std::conj(matrix[2]) + matrix[1] * std::conj(matrix[3]),
		matrix[2] * std::conj(matrix[0]) + matrix[3] * std::conj(matrix[1]),
		matrix[2] * std::conj(matrix[2]) + matrix[3] * std::conj(matrix[3])
	};
	// Use quadratic formula, with a=1.
       double
	 b = -temp[0].real() - temp[3].real(),
	 c = temp[0].real()*temp[3].real() - (temp[1]*temp[2]).real(),
	 d = b*b - (4.0*1.0)*c;
	double
	  sqrtd = sqrt(d);

	e1 = sqrt((-b + sqrtd) * 0.5);
	e2 = sqrt((-b - sqrtd) * 0.5);
}
