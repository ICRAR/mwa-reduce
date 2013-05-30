#include "calibrationmethod.h"

CalibrationMethod::CalibrationMethod(size_t nChannels, size_t nAntenna) :
	_data(nChannels, nAntenna),
	_weights(nChannels, nAntenna),
	_nChannels(nChannels),
	_nAntenna(nAntenna)
{
}

void CalibrationMethod::AddData(const std::complex<float>* data, const float* weights, const std::complex<double>* predictedValues, size_t antenna1, size_t antenna2)
{
	std::complex<double> *destDataPtr = _data.ValuePtr(antenna1, antenna2);
	double *destWeightPtr = _weights.ValuePtr(antenna1, antenna2);
	
	const std::complex<float>* dataEndPtr = data + _nChannels * 4;
	
	std::complex<double> model[4];
	while(data != dataEndPtr)
	{
		for(size_t i=0; i!=4; ++i)
		{
			*destDataPtr = std::complex<double>(*data);
			model[i] = *predictedValues;
			*destWeightPtr = *weights;
			
			++data;
			++weights;
			++predictedValues;
			++destDataPtr;
			++destWeightPtr;
		}
		
		multiplyWithInverse2x2(destDataPtr-4, model);
	}
}

void CalibrationMethod::multiplyWithInverse2x2(std::complex<double>* lhs, std::complex<double>* rhs)
{
	std::complex<double> oneOverDeterminant = 1.0 / ((rhs[0]*rhs[3]) - (rhs[1]*rhs[2]));
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
