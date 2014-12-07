#include "imageweights.h"
#include "banddata.h"
#include "multibanddata.h"

#include "msprovider/msprovider.h"

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights(const WeightMode& weightMode, size_t imageWidth, size_t imageHeight, double pixelScaleX, double pixelScaleY, double superWeight) :
	_weightMode(weightMode),
	_imageWidth(round(double(imageWidth) / superWeight)),
	_imageHeight(round(double(imageHeight) / superWeight)),
	_pixelScaleX(pixelScaleX),
	_pixelScaleY(pixelScaleY),
	_totalSum(0.0),
	_isGriddingFinished(false)
{
	if(_imageWidth%2 != 0) ++_imageWidth;
	if(_imageHeight%2 != 0) ++_imageHeight;
	_grid.assign(_imageWidth*_imageHeight/2, 0.0);
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

void ImageWeights::Grid(casa::MeasurementSet& ms, const MSSelection& selection)
{
	if(_isGriddingFinished)
		throw std::runtime_error("Grid() called after a call to FinishGridding()");
	const MultiBandData bandData(ms.spectralWindow(), ms.dataDescription());
	casa::ROScalarColumn<int> antenna1Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> antenna2Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROScalarColumn<int> fieldIdColumn(ms, casa::MS::columnName(casa::MSMainEnums::FIELD_ID));
	casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
	casa::ROArrayColumn<double> uvwColumn(ms, casa::MS::columnName(casa::MSMainEnums::UVW));
	casa::ROArrayColumn<float> weightColumn(ms, casa::MS::columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, casa::MS::columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casa::MSMainEnums::DATA_DESC_ID));
	
	const casa::IPosition shape(flagColumn.shape(0));
	
	bool isWeightDefined = weightColumn.isDefined(0);
	bool hasWeights = false;
	if(isWeightDefined)
	{
		casa::IPosition modelShape = weightColumn.shape(0);
		hasWeights = (modelShape == shape);
	}
		
	const size_t polarizationCount = shape[0];
	
	casa::Array<casa::Complex> dataArr(shape);
	casa::Array<bool> flagArr(shape);
	casa::Array<float> weightArr(shape);
	size_t timestep = 0;
	double time = timeColumn(0);
	
	if(!hasWeights)
		weightArr.set(1.0);
	
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		const int a1 = antenna1Column(row), a2 = antenna2Column(row), fieldId = fieldIdColumn(row);
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		const casa::Vector<double> uvw = uvwColumn(row);
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvw))
		{
			flagColumn.get(row, flagArr);
			if(hasWeights)
				weightColumn.get(row, weightArr);
			const BandData& curBand = bandData[dataDescIdColumn(row)];
			
			bool* flagIter = flagArr.cbegin();
			float* weightIter = weightArr.cbegin();
			
			double uInM = uvw(0), vInM = uvw(1);
			if(vInM < 0.0)
			{
				uInM = -uInM;
				vInM = -vInM;
			}
			
			size_t startChannel, endChannel;
			if(selection.HasChannelRange())
			{
				startChannel = selection.ChannelRangeStart();
				endChannel = selection.ChannelRangeEnd();
			}
			else {
				startChannel = 0;
				endChannel = curBand.ChannelCount();
			}
	
			for(size_t ch=startChannel; ch!=endChannel; ++ch)
			{
				double
					u = uInM / curBand.ChannelWavelength(ch),
					v = vInM / curBand.ChannelWavelength(ch);
				double x = round(u*_imageWidth*_pixelScaleX + _imageWidth/2);
				double y = round(v*_imageHeight*_pixelScaleY);
					
				if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
				{
					for(size_t p=0; p!=polarizationCount; ++p)
					{
						if(!*flagIter)
						{
								size_t index = (size_t) x + (size_t) y*_imageWidth;
								_grid[index] += *weightIter;
								_totalSum += *weightIter;
						}
						++flagIter;
						++weightIter;
					}
				}
				else {
					for(size_t p=0; p!=polarizationCount; ++p)
					{
						++flagIter;
						++weightIter;
					}
				}
			}
		}
	}
}

void ImageWeights::Grid(MSProvider& msProvider, const MSSelection& selection)
{
	if(_isGriddingFinished)
		throw std::runtime_error("Grid() called after a call to FinishGridding()");
	if(_weightMode.RequiresGridding())
	{
		const MultiBandData bandData(msProvider.MS().spectralWindow(), msProvider.MS().dataDescription());
		MultiBandData selectedBand;
		if(selection.HasChannelRange())
			selectedBand = MultiBandData(bandData, selection.ChannelRangeStart(), selection.ChannelRangeEnd());
		else
			selectedBand = bandData;
		std::vector<float> weightBuffer(selectedBand.MaxChannels());
		
		msProvider.Reset();
		do
		{
			double uInM, vInM, wInM;
			size_t dataDescId;
			msProvider.ReadMeta(uInM, vInM, wInM, dataDescId);
			msProvider.ReadWeights(weightBuffer.data());
			const BandData& curBand = selectedBand[dataDescId];
			if(vInM < 0.0)
			{
				uInM = -uInM;
				vInM = -vInM;
			}
			
			const float* weightIter = weightBuffer.data();
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				double
					u = uInM / curBand.ChannelWavelength(ch),
					v = vInM / curBand.ChannelWavelength(ch);
				double x = round(u*_imageWidth*_pixelScaleX + _imageWidth/2);
				double y = round(v*_imageHeight*_pixelScaleY);
					
				if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
				{
					size_t index = (size_t) x + (size_t) y*_imageWidth;
					_grid[index] += *weightIter;
					_totalSum += *weightIter;
				}
				++weightIter;
			}
		} while(msProvider.NextRow());
	}
}

void ImageWeights::FinishGridding()
{
	if(_isGriddingFinished)
		throw std::runtime_error("FinishGridding() called twice");
	_isGriddingFinished = true;
	
	switch(_weightMode.Mode())
	{
		case WeightMode::BriggsWeighted:
		{
			double avgW = 0.0;
			for(ao::uvector<double>::const_iterator i=_grid.begin(); i!=_grid.end(); ++i)
				avgW += *i * *i;
			avgW /= _totalSum;
			double numeratorSqrt = 5.0 * exp10(-_weightMode.BriggsRobustness());
			double sSq = numeratorSqrt*numeratorSqrt / avgW;
			for(ao::uvector<double>::iterator i=_grid.begin(); i!=_grid.end(); ++i)
			{
				*i = 1.0 / (1.0 + *i * sSq);
			}
		}
		break;
		case WeightMode::UniformWeighted:
		{
			for(ao::uvector<double>::iterator i=_grid.begin(); i!=_grid.end(); ++i)
			{
				if(*i != 0.0)
					*i = 1.0 / *i;
				else
					*i = 0.0;
			}
		}
		break;
		case WeightMode::NaturalWeighted:
		{
			_grid.assign(_imageWidth*_imageHeight/2, 1.0);
		}
		break;
		case WeightMode::DistanceWeighted:
		{
			ao::uvector<double>::iterator i = _grid.begin();
			for(size_t y=0; y!=_imageHeight/2; ++y)
			{
				for(size_t x=0; x!=_imageWidth; ++x)
				{
					double u = double(x-_imageWidth/2) / (_imageWidth*_pixelScaleX);
					double v = double(y) / (_imageHeight*_pixelScaleY);
					*i = GetInverseTaperedWeight(u, v);
					++i;
				}
			}
		}
		break;
	}
}

void ImageWeights::Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep)
{
	if(_isGriddingFinished)
		throw std::runtime_error("Grid() called after a call to FinishGridding()");
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
			double x = round(uTimesLambda*_imageWidth*_pixelScaleX/wavelength + _imageWidth/2);
			double y = round(vTimesLambda*_imageHeight*_pixelScaleY/wavelength);
			if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
			{
				size_t index = (size_t) x + (size_t) y*_imageWidth;
				_grid[index] += 1.0;
			}
		}
	}
}

void ImageWeights::SetMinUVRange(double minUVInLambda)
{
	ao::uvector<double>::iterator i = _grid.begin();
	const double minSq = minUVInLambda*minUVInLambda;
	int halfWidth = _imageWidth/2;
	for(size_t y=0; y!=_imageHeight/2; ++y)
	{
		for(size_t x=0; x!=_imageWidth; ++x)
		{
			int xi = int(x)-halfWidth;
			double u = double(xi) / (_imageWidth*_pixelScaleX);
			double v = double(y) / (_imageHeight*_pixelScaleY);
			if(u*u + v*v < minSq)
				*i = 0.0;
			++i;
		}
	}
}

void ImageWeights::SetMaxUVRange(double maxUVInLambda)
{
	ao::uvector<double>::iterator i = _grid.begin();
	const double maxSq = maxUVInLambda*maxUVInLambda;
	int halfWidth = _imageWidth/2;
	for(size_t y=0; y!=_imageHeight/2; ++y)
	{
		for(size_t x=0; x!=_imageWidth; ++x)
		{
			int xi = int(x)-halfWidth;
			double u = double(xi) / (_imageWidth*_pixelScaleX);
			double v = double(y) / (_imageHeight*_pixelScaleY);
			if(u*u + v*v > maxSq)
				*i = 0.0;
			++i;
		}
	}
}
