#include "fitsreader.h"
#include "fitswriter.h"
#include "imagecoordinates.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <string.h>

class ImageInfo
{
	public:
		ImageInfo() :
			width(0), height(0),
			ra(0.0), dec(0.0),
			pixelSizeX(0.0), pixelSizeY(0.0)
		{
		}
		ImageInfo(size_t size) :
			values(size), weights(size),
			width(0), height(0),
			ra(0.0), dec(0.0),
			pixelSizeX(0.0), pixelSizeY(0.0)
		{
		}
		
		std::vector<double> values, weights;
		size_t width, height;
		double ra, dec, pixelSizeX, pixelSizeY;
};

void Regrid(ImageInfo &destImage, const ImageInfo &sourceImage)
{
	double *outImagePtr = &destImage.values[0], *outWeightPtr = &destImage.weights[0];
	size_t withinField = 0;
	for(size_t y=0; y!=destImage.height; ++y)
	{
		for(size_t x=0; x!=destImage.width; ++x)
		{
			long double l, m;
			ImageCoordinates::FitsXYToLM<long double>(x, y, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, l, m);
			long double ra, dec;
			ImageCoordinates::LMToRaDec<long double>(l, m, destImage.ra, destImage.dec, ra, dec);
			
			long double sourceL, sourceM;
			ssize_t sourceX, sourceY;
			ImageCoordinates::RaDecToLM<long double>(ra, dec, sourceImage.ra, sourceImage.dec, sourceL, sourceM);
			ImageCoordinates::FitsLMToXY<long double>(sourceL, sourceM, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, sourceX, sourceY);
			
			if(sourceX >= 0 && sourceX < (ssize_t) sourceImage.width && sourceY >= 0 && sourceY < (ssize_t) sourceImage.height)
			{
				size_t sourceValueIndex = sourceY * sourceImage.width + sourceX;
				double sourceValue = sourceImage.values[sourceValueIndex];
				double sourceWeight = sourceImage.weights[sourceValueIndex];
				
				*outImagePtr += sourceValue * sourceWeight;
				*outWeightPtr += sourceWeight;
				
				++withinField;
			}
			
			++outImagePtr;
			++outWeightPtr;
		}
	}
	std::cout << (withinField * 100 / (destImage.height*destImage.width)) << "% pixels fitted in new image.\n";
}

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cerr << "Syntax: regridimg [options] <outimage> <outweights> <template> <inpimage1> <inpweights1> [<inpimage2> <inpweights2> ...]\n"
			"All images should be fits files. First image will define the size and pointing centre.\n";
	}
	int argi = 1;
	bool overrideSize = false;
	size_t width, height;
	if(strcmp(argv[argi], "-s") == 0)
	{
		overrideSize = true;
		width = atoi(argv[argi+1]);
		height = atoi(argv[argi+2]);
		argi += 3;
	}
	const char
		*outImageName = argv[argi],
		*outWeightName = argv[argi+1],
		*templateName = argv[argi+2];
	
	FitsReader templateReader(templateName);
	if(!overrideSize)
	{
		width = templateReader.ImageWidth();
		height = templateReader.ImageHeight();
	}
	const size_t size = width * height;
	ImageInfo outImage(size);
	outImage.width = width;
	outImage.height = height;
	outImage.ra = templateReader.PhaseCentreRA();
	outImage.dec = templateReader.PhaseCentreDec();
	outImage.pixelSizeX = templateReader.PixelSizeX();
	outImage.pixelSizeY = templateReader.PixelSizeY();
	for(size_t i=0; i!=size; ++i)
	{
		outImage.values[i] = 0.0;
		outImage.weights[i] = 0.0;
	}
	argi += 3;
	for(; argi + 1 < argc; argi += 2)
	{
		const char *inpImageName = argv[argi];
		const char *inpWeightName = argv[argi+1];
		
		std::cout << "Regridding " << inpImageName << "...\n";
		
		FitsReader inpReader(inpImageName);
		FitsReader weightsReader(inpWeightName);
		if(weightsReader.ImageWidth() != inpReader.ImageWidth() || weightsReader.ImageHeight() != inpReader.ImageHeight())
			throw std::runtime_error("Weights and image do not have same size");
		
		ImageInfo inpImage(inpReader.ImageWidth() * inpReader.ImageHeight());
		
		inpReader.Read<double>(&inpImage.values[0]);
		weightsReader.Read<double>(&inpImage.weights[0]);
		
		inpImage.width = inpReader.ImageWidth(),
		inpImage.height = inpReader.ImageHeight();
		inpImage.ra = inpReader.PhaseCentreRA();
		inpImage.dec = inpReader.PhaseCentreDec();
		inpImage.pixelSizeX = inpReader.PixelSizeX();
		inpImage.pixelSizeY = inpReader.PixelSizeY();
		
		Regrid(outImage, inpImage);
	}
	
	// Divide the weight out
	std::cout << "Applying weights...\n";
	const double *weightsIter = &outImage.weights[0];
	double *imageEnd = &outImage.values[0] + (outImage.width * outImage.height);
	for(double *imagePtr=&outImage.values[0]; imagePtr!=imageEnd; ++imagePtr)
	{
		*imagePtr /=  *weightsIter;
		++weightsIter;
	}
	
	std::cout << "Writing " << outImageName << "...\n";
	FitsWriter imgWriter(outImageName);
	imgWriter.Write<double>(&outImage.values[0], outImage.width, outImage.height, outImage.ra, outImage.dec, outImage.pixelSizeX, outImage.pixelSizeY);
	
	std::cout << "Writing " << outWeightName << "...\n";
	FitsWriter weightsWriter(outWeightName);
	weightsWriter.Write<double>(&outImage.weights[0], outImage.width, outImage.height, outImage.ra, outImage.dec, outImage.pixelSizeX, outImage.pixelSizeY);
}
