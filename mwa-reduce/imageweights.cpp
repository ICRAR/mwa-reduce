#include "imageweights.h"
#include "banddata.h"
#include "multibanddata.h"

#include "msprovider/msprovider.h"

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights(size_t imageWidth, size_t imageHeight, double pixelScaleX, double pixelScaleY, double superWeight) :
	_imageWidth(round(double(imageWidth) / superWeight)),
	_imageHeight(round(double(imageHeight) / superWeight)),
	_pixelScaleX(pixelScaleX),
	_pixelScaleY(pixelScaleY)
{
	if(_imageWidth%2 != 0) ++_imageWidth;
	if(_imageHeight%2 != 0) ++_imageHeight;
	_sum.assign(_imageWidth*_imageHeight/2, 0.0);
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

void ImageWeights::Grid(casa::MeasurementSet& ms, WeightMode weightMode, const MSSelection& selection)
{
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
	double totalSum = 0.0;
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
								_sum[index] += *weightIter;
								if(weightMode.IsBriggs())
									totalSum += *weightIter;
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
	// TODO this is not right: can't update Sum here since there might be more mses be gridded
	if(weightMode.IsBriggs())
	{
		double avgW = 0.0;
		for(ao::uvector<double>::const_iterator i=_sum.begin(); i!=_sum.end(); ++i)
			avgW += *i * *i;
		avgW /= totalSum;
		double numeratorSqrt = 5.0 * exp10(-weightMode.BriggsRobustness());
		double sSq = numeratorSqrt*numeratorSqrt / avgW;
		for(ao::uvector<double>::iterator i=_sum.begin(); i!=_sum.end(); ++i)
		{
			*i = 1.0 / (1.0 + *i * sSq);
		}
	}
}

void ImageWeights::Grid(MSProvider& msProvider, WeightMode weightMode, const MSSelection& selection)
{
	const MultiBandData bandData(msProvider.MS().spectralWindow(), msProvider.MS().dataDescription());
	MultiBandData selectedBand;
	if(selection.HasChannelRange())
		selectedBand = MultiBandData(bandData, selection.ChannelRangeStart(), selection.ChannelRangeEnd());
	else
		selectedBand = bandData;
	double totalSum = 0.0;
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
				_sum[index] += *weightIter;
				if(weightMode.IsBriggs())
					totalSum += *weightIter;
			}
			++weightIter;
		}
	} while(msProvider.NextRow());
	
	// TODO this is not right: can't update Sum here since there might be more mses be gridded
	if(weightMode.IsBriggs())
	{
		double avgW = 0.0;
		for(ao::uvector<double>::const_iterator i=_sum.begin(); i!=_sum.end(); ++i)
			avgW += *i * *i;
		avgW /= totalSum;
		double numeratorSqrt = 5.0 * exp10(-weightMode.BriggsRobustness());
		double sSq = numeratorSqrt*numeratorSqrt / avgW;
		for(ao::uvector<double>::iterator i=_sum.begin(); i!=_sum.end(); ++i)
		{
			*i = 1.0 / (1.0 + *i * sSq);
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
			double x = round(uTimesLambda*_imageWidth*_pixelScaleX/wavelength + _imageWidth/2);
			double y = round(vTimesLambda*_imageHeight*_pixelScaleY/wavelength);
			if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
			{
				size_t index = (size_t) x + (size_t) y*_imageWidth;
				_sum[index] += 1.0;
			}
		}
	}
}
