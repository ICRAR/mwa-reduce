#include "fitsreader.h"
#include "fitswriter.h"
#include "imagecoordinates.h"
#include "radeccoord.h"

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

void getBoundingBox(const ImageInfo &destImage, const ImageInfo &sourceImage, size_t &destXLeft, size_t &destYTop, size_t &destXRight, size_t &destYBottom)
{
	long double l, m;
	long double ra, dec;
	long double destL, destM;
	int destX, destY;
	
	// Traverse top and bottom edges and find outer coordinates for those
	for(size_t x=0; x!=sourceImage.width; ++x)
	{
		// Determine RA,DEC for top edge position in source image
		ImageCoordinates::XYToLM<long double>(x, 0, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, l, m);
		ImageCoordinates::LMToRaDec<long double>(l, m, sourceImage.ra, sourceImage.dec, ra, dec);
		// Conv RA,DEC to dest image positions
		ImageCoordinates::RaDecToLM<long double>(ra, dec, destImage.ra, destImage.dec, destL, destM);
		ImageCoordinates::LMToXY<long double>(destL, destM, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, destX, destY);
		
		if(destX >= (int) destImage.width) destX = destImage.width-1;
		if(destX < 0) destX = 0;
		if(destY >= (int) destImage.height) destY = destImage.height-1;
		if(destY < 0) destY = 0;
		
		// Initialize coordinates during first iteration
		if(x == 0)
		{
			destXLeft = destX;
			destXRight = destX;
			destYTop = destY;
			destYBottom = destY;
		}
		
		destXLeft = std::min((size_t) destX, destXLeft);
		destYTop = std::min((size_t) destY, destYTop);
		destXRight = std::max((size_t) destX, destXRight);
		destYBottom = std::max((size_t) destY, destYBottom);
		
		// Determine RA,DEC for bottom edge position in source image
		ImageCoordinates::XYToLM<long double>(x, sourceImage.height-1, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, l, m);
		ImageCoordinates::LMToRaDec<long double>(l, m, sourceImage.ra, sourceImage.dec, ra, dec);
		// Conv RA,DEC to dest image positions
		ImageCoordinates::RaDecToLM<long double>(ra, dec, destImage.ra, destImage.dec, destL, destM);
		ImageCoordinates::LMToXY<long double>(destL, destM, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, destX, destY);
		
		if(destX >= (int) destImage.width) destX = destImage.width-1;
		if(destX < 0) destX = 0;
		if(destY >= (int) destImage.height) destY = destImage.height-1;
		if(destY < 0) destY = 0;
		
		destXLeft = std::min((size_t) destX, destXLeft);
		destYTop = std::min((size_t) destY, destYTop);
		destXRight = std::max((size_t) destX, destXRight);
		destYBottom = std::max((size_t) destY, destYBottom);
	}
	
	// Traverse left and right edges and find outer coordinates for those
	for(size_t y=0; y!=sourceImage.height; ++y)
	{
		// Determine RA,DEC for left edge position in source image
		ImageCoordinates::XYToLM<long double>(0, y, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, l, m);
		ImageCoordinates::LMToRaDec<long double>(l, m, sourceImage.ra, sourceImage.dec, ra, dec);
		// Conv RA,DEC to dest image positions
		ImageCoordinates::RaDecToLM<long double>(ra, dec, destImage.ra, destImage.dec, destL, destM);
		ImageCoordinates::LMToXY<long double>(destL, destM, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, destX, destY);
		
		if(destX >= (int) destImage.width) destX = destImage.width-1;
		if(destX < 0) destX = 0;
		if(destY >= (int) destImage.height) destY = destImage.height-1;
		if(destY < 0) destY = 0;
		
		destXLeft = std::min((size_t) destX, destXLeft);
		destYTop = std::min((size_t) destY, destYTop);
		destXRight = std::max((size_t) destX, destXRight);
		destYBottom = std::max((size_t) destY, destYBottom);
		
		// Determine RA,DEC for right edge position in source image
		ImageCoordinates::XYToLM<long double>(sourceImage.width-1, y, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, l, m);
		ImageCoordinates::LMToRaDec<long double>(l, m, sourceImage.ra, sourceImage.dec, ra, dec);
		// Conv RA,DEC to dest image positions
		ImageCoordinates::RaDecToLM<long double>(ra, dec, destImage.ra, destImage.dec, destL, destM);
		ImageCoordinates::LMToXY<long double>(destL, destM, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, destX, destY);
		
		if(destX >= (int) destImage.width) destX = destImage.width-1;
		if(destX < 0) destX = 0;
		if(destY >= (int) destImage.height) destY = destImage.height-1;
		if(destY < 0) destY = 0;
		
		destXLeft = std::min((size_t) destX, destXLeft);
		destYTop = std::min((size_t) destY, destYTop);
		destXRight = std::max((size_t) destX, destXRight);
		destYBottom = std::max((size_t) destY, destYBottom);
	}
}

void Regrid(ImageInfo &destImage, const ImageInfo &sourceImage)
{
	size_t withinField = 0;
	
	size_t xLeft, xRight, yTop, yBottom;
	getBoundingBox(destImage, sourceImage, xLeft, yTop, xRight, yBottom);
	
	for(size_t y=yTop; y<yBottom; ++y)
	{
		double
			*outImagePtr = &destImage.values[y * destImage.width + xLeft],
			*outWeightPtr = &destImage.weights[y * destImage.width + xLeft];
		for(size_t x=xLeft; x<xRight; ++x)
		{
			long double l, m;
			ImageCoordinates::XYToLM<long double>(x, y, destImage.pixelSizeX, destImage.pixelSizeY, destImage.width, destImage.height, l, m);
			long double ra, dec;
			ImageCoordinates::LMToRaDec<long double>(l, m, destImage.ra, destImage.dec, ra, dec);
			
			long double sourceL, sourceM;
			int sourceX, sourceY;
			ImageCoordinates::RaDecToLM<long double>(ra, dec, sourceImage.ra, sourceImage.dec, sourceL, sourceM);
			ImageCoordinates::LMToXY<long double>(sourceL, sourceM, sourceImage.pixelSizeX, sourceImage.pixelSizeY, sourceImage.width, sourceImage.height, sourceX, sourceY);
			
			if(sourceX >= 0 && sourceX < (int) sourceImage.width && sourceY >= 0 && sourceY < (int) sourceImage.height)
			{
				size_t sourceValueIndex = sourceY * sourceImage.width + sourceX;
				double sourceValue = sourceImage.values[sourceValueIndex];
				double sourceBeam = sourceImage.weights[sourceValueIndex];// * sourceImage.weights[sourceValueIndex];
				
				*outImagePtr += sourceValue * sourceBeam;
				*outWeightPtr += sourceBeam * sourceBeam;
				
				++withinField;
			}
			
			++outImagePtr;
			++outWeightPtr;
		}
	}
	std::cout << "Bounding box: (" << xLeft << ',' << yTop << ")-(" << xRight << ',' << yBottom << "), "
	<< (withinField * 100 / (destImage.height*destImage.width)) << "% pixels fitted in new image.\n";
}

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cerr << "Syntax: regridimg [options] <outimage> <outweights> <template> <inpimage1> <inpweights1> [<inpimage2> <inpweights2> ...]\n"
			"All images should be fits files. First image will define the size and pointing centre.\n";
	}
	int argi = 1;
	bool overrideSize = false, overrideCentre = false;
	size_t width = 0, height = 0;
	long double altRA = 0.0, altDec = 0.0;
	while(argv[argi][0] == '-')
	{
		if(strcmp(argv[argi], "-s") == 0)
		{
			overrideSize = true;
			width = atoi(argv[argi+1]);
			height = atoi(argv[argi+2]);
			argi += 3;
		}
		else if(strcmp(argv[argi], "-c") == 0)
		{
			overrideCentre = true;
			altRA = RaDecCoord::ParseRA(argv[argi+1]);
			altDec = RaDecCoord::ParseDec(argv[argi+2]);
			argi += 3;
		}
		else throw std::runtime_error("Bad parameter");
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
	if(overrideCentre)
	{
		outImage.ra = altRA;
		outImage.dec = altDec;
	} else {
		outImage.ra = templateReader.PhaseCentreRA();
		outImage.dec = templateReader.PhaseCentreDec();
	}
	outImage.pixelSizeX = templateReader.PixelSizeX();
	outImage.pixelSizeY = templateReader.PixelSizeY();
	for(size_t i=0; i!=size; ++i)
	{
		outImage.values[i] = 0.0;
		outImage.weights[i] = 0.0;
	}
	argi += 3;
	int argStart = argi;
	for(; argi + 1 < argc; argi += 2)
	{
		const char *inpImageName = argv[argi];
		const char *inpWeightName = argv[argi+1];
		
		std::cout << "Regridding " << inpImageName << "... (" << ((argi-argStart)/2) << '/' << ((argc-argStart)/2) << ")\n";
		
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
		if((*weightsIter) != 0.0)
			*imagePtr /= *weightsIter;
		else
			*imagePtr = 0.0;
		++weightsIter;
	}
	
	std::cout << "Writing " << outImageName << "...\n";
	FitsWriter imgWriter(outImageName);
	imgWriter.Write<double>(&outImage.values[0], outImage.width, outImage.height, outImage.ra, outImage.dec, outImage.pixelSizeX, outImage.pixelSizeY);
	
	std::cout << "Writing " << outWeightName << "...\n";
	FitsWriter weightsWriter(outWeightName);
	weightsWriter.Write<double>(&outImage.weights[0], outImage.width, outImage.height, outImage.ra, outImage.dec, outImage.pixelSizeX, outImage.pixelSizeY);
}
