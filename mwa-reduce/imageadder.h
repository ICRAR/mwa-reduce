#ifndef IMAGE_ADDER_H
#define IMAGE_ADDER_H

#include "fitsreader.h"
#include "fitswriter.h"
#include "uvector.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <memory>

class ImageAdder
{
public:
	ImageAdder() :
		_width(0), _height(0), _count(0),
		_outImage(0), _outWeights(0),
		_imgWriter(),
		_frequencySum(0.0), _lowestFreq(0.0), _highestFreq(0.0),
		_normalizationFactorSum(0.0), _normalizationFactorWeightSum(0.0),
		_totalImageWeight(0.0), _totalImageWeightWeightSum(0.0)
	{
	}
	
	void Add(const std::string& imageFilename, const std::string& weightFilename="-1", const std::string& imaginaryWeightFilename="")
	{
		FitsReader inpReader(imageFilename);
		if(_outImage.empty())
		{
			_width = inpReader.ImageWidth(),
			_height = inpReader.ImageHeight();
			_imgWriter.reset(new FitsWriter(inpReader));
			
			const size_t size = _width * _height;
			_outImage.assign(size, 0.0);
			_outWeights.assign(size, 0.0);
		} else {
			if(_width != inpReader.ImageWidth() || _height != inpReader.ImageHeight())
				throw std::runtime_error("Not all images have same size");
		}
		
		double thisFrequency = inpReader.Frequency();
		_frequencySum += thisFrequency;
		if(_count == 0)
		{
			_lowestFreq = thisFrequency - inpReader.Bandwidth()*0.5;
			_highestFreq = thisFrequency + inpReader.Bandwidth()*0.5;
		}
		else {
			if(thisFrequency - inpReader.Bandwidth()*0.5 < _lowestFreq)
				_lowestFreq = thisFrequency - inpReader.Bandwidth()*0.5;
			if(thisFrequency + inpReader.Bandwidth()*0.5 > _highestFreq)
				_highestFreq = thisFrequency + inpReader.Bandwidth()*0.5;
		}
		_count++;
		
		ao::uvector<double> inpImage(_width*_height), weightImage(_width*_height);
		
		inpReader.Read<double>(&inpImage[0]);
		double wscImageWeight;
		if(!inpReader.ReadDoubleKeyIfExists("WSCIMGWG", wscImageWeight))
		{
			std::cerr << "Warning: Keyword WSCIMGWG not found!\n";
			wscImageWeight = 1.0;
		}
		double wscNormalizationFactor;
		if(!inpReader.ReadDoubleKeyIfExists("WSCNORMF", wscNormalizationFactor))
		{
			std::cerr << "Warning: Keyword WSCNORM not found!\n";
			wscNormalizationFactor = 1.0;
		}
		
		if(std::string(weightFilename) == "-1")
		{
			weightImage.assign(_width*_height, wscImageWeight);
		}
		else if(!imaginaryWeightFilename.empty())
		{
			FitsReader realReader(weightFilename), imagReader(imaginaryWeightFilename);
			if(realReader.ImageWidth() != _width || realReader.ImageHeight() != _height)
				throw std::runtime_error("Real beam and image do not have same size");
			if(imagReader.ImageWidth() != _width || imagReader.ImageHeight() != _height)
				throw std::runtime_error("Imaginary beam and image do not have same size");
			
			ao::uvector<double> realImage(_width*_height), imagImage(_width*_height);
			realReader.Read<double>(&realImage[0]);
			imagReader.Read<double>(&imagImage[0]);
			for(size_t j=0; j!=_width*_height; ++j)
			{
				double r = realImage[j], i = imagImage[j];
				weightImage[j] = (r*r + i*i) * wscImageWeight;
			}
		}
		else {
			FitsReader weightsReader(weightFilename);
			if(weightsReader.ImageWidth() != _width || weightsReader.ImageHeight() != _height)
				throw std::runtime_error("Weights and image do not have same size");
			weightsReader.Read<double>(&weightImage[0]);
			for(size_t j=0; j!=_width*_height; ++j)
				weightImage[j] *= wscImageWeight;
		}
		
		double centralWeight = weightImage[_width/2 + (_height/2)*_width];
		_normalizationFactorSum += wscNormalizationFactor * wscImageWeight * centralWeight;
		_normalizationFactorWeightSum += wscImageWeight * centralWeight;
		_totalImageWeight += wscImageWeight * centralWeight;
		_totalImageWeightWeightSum += centralWeight;
		
		// Add the images in
		double *outImagePtr = _outImage.data(), *outWeightPtr = _outWeights.data();
		ao::uvector<double>::iterator inpWeightsIter = weightImage.begin();
		for(ao::uvector<double>::iterator i=inpImage.begin(); i!=inpImage.end(); ++i)
		{
			double beamVal = *inpWeightsIter;
			*outImagePtr +=  (*i) * beamVal;
			*outWeightPtr += beamVal * beamVal;
			
			++inpWeightsIter;
			++outImagePtr;
			++outWeightPtr;
		}
		
		std::cout << '.' << std::flush;
	}
	
	void Finish(const std::string& outputFilename, const std::string& outputWeightsFilename="")
	{
		// Divide the weight out
		const double *weightsIter = _outWeights.data();
		double *imageEnd = _outImage.data() + (_width * _height);
		for(double *imagePtr=_outImage.data(); imagePtr!=imageEnd; ++imagePtr)
		{
			double weight = *weightsIter;
			*imagePtr /=  weight;
			
			++weightsIter;
		}
		
		_imgWriter->SetFrequency(_frequencySum / _count, (_highestFreq - _lowestFreq));
		_imgWriter->SetExtraKeyword("WSCNORMF", _normalizationFactorSum / _normalizationFactorWeightSum);
		_imgWriter->SetExtraKeyword("WSCIMGWG", _totalImageWeight / _totalImageWeightWeightSum);
		_imgWriter->Write<double>(outputFilename, _outImage.data());
		
		if(!outputWeightsFilename.empty() && outputWeightsFilename!="-")
			_imgWriter->Write<double>(outputWeightsFilename, _outWeights.data());
	}
private:
	size_t _width, _height, _count;
	ao::uvector<double> _outImage, _outWeights;
	std::unique_ptr<FitsWriter> _imgWriter;
	double _frequencySum, _lowestFreq, _highestFreq;
	double _normalizationFactorSum, _normalizationFactorWeightSum;
	double _totalImageWeight, _totalImageWeightWeightSum;
};

#endif
