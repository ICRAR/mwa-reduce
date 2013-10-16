#include "imageweights.h"
#include "banddata.h"
#include "multibanddata.h"

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights(size_t imageWidth, size_t imageHeight, double pixelScale) : 
	_imageWidth(imageWidth),
	_imageHeight(imageHeight),
	_pixelScale(pixelScale),
	_sum(_imageWidth*_imageHeight/2)
{
}

double ImageWeights::GetUniformWeight(double u, double v)
{
	if(v < 0.0) {
		u = -u;
		v = -v;
	}
	double x = round(u*_imageWidth*_pixelScale + _imageWidth/2);
	double y = round(v*_imageHeight*_pixelScale);
	if(x >= 0.0 && x < _imageWidth && y < _imageHeight)
		return 1.0 / (double) _sum[(size_t) x + (size_t) y*_imageWidth];
	else {
		return 0.0;
	}
}

double ImageWeights::ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep)
{
	double weightSum = 0.0;
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		if(flags[ch])
		{
			data[ch] = 0.0;
		} else
		{
			double wavelength = frequencyToWavelength(lowestFrequency + frequencyStep*ch);
			double u = uTimesLambda/wavelength;
			double v = vTimesLambda/wavelength;
			double weight = GetWeight(u, v);
			weightSum += weight;
			data[ch] *= weight;
		}
	}
	return weightSum / channelCount;
}

void ImageWeights::Grid(casa::MeasurementSet& ms)
{
	const MultiBandData bandData(ms.spectralWindow(), ms.dataDescription());
	casa::ROScalarColumn<int> antenna1Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> antenna2Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, casa::MS::columnName(casa::MSMainEnums::UVW));
	casa::ROArrayColumn<float> weightColumn(ms, casa::MS::columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, casa::MS::columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casa::MSMainEnums::DATA_DESC_ID));
	
	const casa::IPosition shape(flagColumn.shape(0));
	const size_t polarizationCount = shape[0];
	
	casa::Array<casa::Complex> dataArr(shape);
	casa::Array<bool> flagArr(shape);
	casa::Array<float> weightArr(shape);
	
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		const int a1 = antenna1Column(row), a2 = antenna2Column(row);
		if(a1 != a2)
		{
			flagColumn.get(row, flagArr);
			weightColumn.get(row, weightArr);
			const casa::Vector<double> uvw = uvwColumn(row);
			const BandData& curBand = bandData[dataDescIdColumn(row)];
			
			bool* flagIter = flagArr.cbegin();
			float* weightIter = weightArr.cbegin();
			
			double uInM = uvw(0), vInM = uvw(1);
			if(vInM < 0.0)
			{
				uInM = -uInM;
				vInM = -vInM;
			}
			
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double
					u = uInM / curBand.ChannelWavelength(ch),
					v = vInM / curBand.ChannelWavelength(ch);
				double x = round(u*_imageWidth*_pixelScale + _imageWidth/2);
				double y = round(v*_imageHeight*_pixelScale);
					
				if(x >= 0.0 && x < _imageWidth && y < _imageHeight)
				{
					for(size_t p=0; p!=polarizationCount; ++p)
					{
						if(!*flagIter)
						{
								size_t index = (size_t) x + (size_t) y*_imageWidth;
								_sum[index] += *weightIter;
						}
						++flagIter;
						++weightIter;
					}
				}
			}
		}
	}
}

void ImageWeights::Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep)
{
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		if(!flags[ch])
		{
			if(vTimesLambda < 0.0)
			{
				uTimesLambda = -uTimesLambda;
				vTimesLambda = -vTimesLambda;
			}
			
			double wavelength = frequencyToWavelength(lowestFrequency + frequencyStep*ch);
			double x = round(uTimesLambda*_imageWidth*_pixelScale/wavelength + _imageWidth/2);
			double y = round(vTimesLambda*_imageHeight*_pixelScale/wavelength);
			if(x >= 0.0 && x < _imageWidth && y < _imageHeight)
			{
				size_t index = (size_t) x + (size_t) y*_imageWidth;
				_sum[index] += 1.0;
			}
		}
	}
}
