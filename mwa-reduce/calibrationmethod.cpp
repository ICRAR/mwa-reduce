#include "calibrationmethod.h"

#include <iostream>
#include <cmath>

CalibrationMethod::CalibrationMethod(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
	_data(nChannels, nAntenna, nTimesteps),
	_model(nChannels, nAntenna, nTimesteps),
	_weights(nChannels, nAntenna, nTimesteps),
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
	for(size_t timestep=0; timestep!=_nTimesteps; ++timestep)
	{
		for(size_t antenna2 = 0; antenna2!=_nAntenna; ++antenna2)
		{
			for(size_t antenna1 = 0; antenna1!=_nAntenna; ++antenna1)
			{
				std::complex<double> *dataPtr = _data.ValuePtr(antenna1, antenna2, timestep);
				double *weightPtr = _weights.ValuePtr(antenna1, antenna2, timestep);
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
					for(size_t p=0; p!=4; ++p)
					{
						*dataPtr *= *weightPtr;
					}
					++weightPtr;
				}
			}
		}
	}
}

void CalibrationMethod::Execute(double precisionLimit)
{
	std::vector<std::complex<double> > nextJones(_nAntenna * 4 * _nChannels);
	
	bool continueIterating;
	
	applyWeightsToData();
	
	do
	{
		for(size_t ant=0; ant!=_nAntenna; ++ant)
		{
			// TODO weight the data
			
			calculateNextIter(ant, &nextJones[ant*4]);
			
			// TODO take weights out
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
				
				if(ant==1 && ch==15)
				{
					std::cout << "Current value of Jones matrix for ant 1, ch 15:\n"
						" (" << jonesPtr[-4] << " " << jonesPtr[-3] << ")\n"
						" (" << jonesPtr[-2] << " " << jonesPtr[-1] << ")\n";
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
		
	} while(continueIterating);
}

void CalibrationMethod::calculateNextIter(size_t ant, std::complex<double> *nextJones)
{
	// Calculate the next Jones estimate for the given ant, based on
	// the previous Jones estimates.
	
	// Calculate:
	// JONES[ant] = SUM_{k!=ant in Nant} MODEL[ant,k] JONES[k] DATA^H[ant,k] times
	//            ( SUM_{k!=ant in Nant} DATA[ant,k] JONES^H[k] JONES[k] DATA^H[ant,k] )^-1
	// (From Mitchel et al., 2008)
	
	for(size_t i=0; i!=4 * _nChannels; ++i)
		nextJones[i] = std::complex<double>(0.0, 0.0);
	
	std::vector<std::complex<double> > rTerm(4 * _nChannels);
	for(size_t t=0; t!=_nTimesteps; ++t)
	{
		for(size_t k=0; k!=_nAntenna; ++k)
		{
			if(k != ant)
			{
				std::complex<double> *dataPtr = _data.ValuePtr(ant, k, t);
				std::complex<double> *modelPtr = _model.ValuePtr(ant, k, t);
				std::complex<double> *jonesPtr = &_jonesSolutions[k * 4];
				std::complex<double> *rTermPtr = &rTerm[0];
				std::complex<double> *nextJonesPtr = nextJones;
			
				for(size_t ch=0; ch!=_nChannels; ++ch)
				{
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
		multiplyWithInverse2x2(&nextJones[ch * 4], &rTerm[ch * 4]);
	}
}

void CalibrationMethod::multiplyWithInverse2x2(std::complex<double>* lhs, std::complex<double>* rhs)
{
	std::complex<double> d = ((rhs[0]*rhs[3]) - (rhs[1]*rhs[2]));
	std::complex<double> oneOverDeterminant = (d != 0.0) ? (1.0 / d) : 0.0;
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
}
