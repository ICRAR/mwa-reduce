#include "imageweights.h"

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights(size_t imageSize, size_t channelCount, double pixelScale, double lowestFrequency, double frequencyStep) : 
	_imageSize(imageSize),
	_channelCount(channelCount),
	_pixelScale(pixelScale),
	_lowestFrequency(lowestFrequency),
	_frequencyStep(frequencyStep)
{
	_counts = new int[_imageSize*_imageSize];
	_sum = new double[_imageSize*_imageSize];
	_sumSq = new double[_imageSize*_imageSize];
	memset(_counts, 0, _imageSize*_imageSize);
}

double ImageWeights::GetCountWeight(double u, double v)
{
	double x = round(u*_imageSize*_pixelScale + _imageSize/2);
	double y = round(v*_imageSize*_pixelScale + _imageSize/2);
	if(x >= 0.0 && x < _imageSize && y >= 0.0 && y < _imageSize)
		return 1.0 / (double) _counts[(size_t) x + (size_t) y*_imageSize];
	else {
		std::cout << "Weights did not fit in grid.\n";
		return 0.0;
	}
}

double ImageWeights::GetUniformWeight(double u, double v)
{
	double x = round(u*_imageSize*_pixelScale + _imageSize/2);
	double y = round(v*_imageSize*_pixelScale + _imageSize/2);
	if(x >= 0.0 && x <= _imageSize && y >= 0.0 && y <= _imageSize)
	{
		size_t index = (size_t) x + (size_t) y*_imageSize;
		double countFact = 1.0 / _counts[index];
		double sum = _sum[index];
		double sumSq = _sumSq[index];
		double sumMeanSquared = sum * sum * countFact;
		double variance = (sumSq - sumMeanSquared) * countFact;
		if(variance!=0.0)
			return 1.0 / variance;
		else
			return 1.0;
	} else {
		std::cout << "Weights did not fit in grid.\n";
		return 0.0;
	}
}

double ImageWeights::ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda)
{
	double weightSum = 0.0;
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		if(flags[ch])
		{
			data[ch] = 0.0;
		} else
		{
			double wavelength = frequencyToWavelength(_lowestFrequency + _frequencyStep*ch);
			double u = uTimesLambda/wavelength;
			double v = vTimesLambda/wavelength;
			double weight = GetWeight(u, v);
			weightSum += weight;
			data[ch] *= weight;
		}
	}
	return weightSum / _channelCount;
}

void ImageWeights::Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda)
{
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		if(!flags[ch])
		{
			double wavelength = frequencyToWavelength(_lowestFrequency + _frequencyStep*ch);
			double x = round(uTimesLambda*_imageSize*_pixelScale/wavelength + _imageSize/2);
			double y = round(vTimesLambda*_imageSize*_pixelScale/wavelength + _imageSize/2);
			if(x >= 0.0 && x < _imageSize && y >= 0.0 && y < _imageSize)
			{
				size_t index = (size_t) x + (size_t) y*_imageSize;
				_counts[index] += 2;
				_sum[index] += data[ch].real() + data[ch].imag();
				_sumSq[index] += data[ch].real()*data[ch].real() + data[ch].imag()*data[ch].imag();
			}
		}
	}
}
