#include "fitsreader.h"
#include "fitswriter.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cmath>

double FindPeak(const std::vector<double> &image, size_t width, size_t &x, size_t &y, bool allowNegativeComponents)
{
	double peakMax = fabs(*image.begin());
	size_t peakIndex = 0;
	size_t index = 0;
	for(std::vector<double>::const_iterator i=image.begin(); i!=image.end(); ++i)
	{
		double value = *i;
		if(allowNegativeComponents) value = fabs(value);
		if(value > peakMax)
		{
			peakIndex = index;
			peakMax = fabs(*i);
		}
		++index;
	}
	x = peakIndex % width;
	y = peakIndex / width;
	return image[x + y*width];
}

void subtractPsf(std::vector<double> &image, const std::vector<double> &psf, size_t width, size_t height, size_t x, size_t y, double factor)
{
	size_t startX, startY, endX, endY;
	ssize_t offsetX = x - width/2, offsetY = y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > 0)
		startY = offsetY;
	else
		startY = 0;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	endY = y + height/2;
	if(endY > height) endY = height;
	
	for(size_t ypos = startY; ypos != endY; ++ypos)
	{
		std::vector<double>::iterator imageIter = image.begin() + ypos * width + startX;
		std::vector<double>::const_iterator psfIter = psf.begin() + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; ++xpos)
		{
			*imageIter -= *psfIter * factor;
			++imageIter;
			++psfIter;
		}
	}
}

int main(int argc, char *argv[])
{
	if(argc != 5)
	{
		std::cerr << "Syntax: imgclean <inpimage> <psfimage> <outpimage> <modelimage>\n"
			"All images should be fits files.\n";
		return -1;
	}
	const char *inpImageName = argv[1];
	const char *psfImageName = argv[2];
	const char *outImageName = argv[3];
	const char *modelImageName = argv[4];
	double threshold = 5.0, subtractionFactor = 0.25;
	size_t maxIter = 50000;
	bool allowNegativeComponents = true, stopWhenDiverging = false;
	
	FitsReader inpReader(inpImageName), psfReader(psfImageName);
	const size_t width = inpReader.ImageWidth();
	const size_t height = inpReader.ImageHeight();
		
	const size_t size = width * height;
	std::vector<double>
		image(size),
		psf(size),
		model(size, 0.0);
	
	inpReader.Read<double>(&image[0]);
	psfReader.Read<double>(&psf[0]);
	
	size_t componentX, componentY;
	double peak = FindPeak(image, width, componentX, componentY, allowNegativeComponents);
	double lastPeak = peak;
	std::cout << "Initial peak: " << peak << '\n';
	size_t iterationNumber = 0;
	while(fabs(peak) > threshold && iterationNumber < maxIter && (!stopWhenDiverging || peak<=lastPeak))
	{
		if(iterationNumber % 10 == 0)
			std::cout << "Iteration " << iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		subtractPsf(image, psf, width, height, componentX, componentY, subtractionFactor * peak);
		model[componentX + componentY*width] += subtractionFactor * peak;
		
		lastPeak = peak;
		peak = FindPeak(image, width, componentX, componentY, allowNegativeComponents);
		++iterationNumber;
	}
	std::cout << "Stopped on peak " << peak << '\n';
	
	FitsWriter imgWriter(outImageName);
	imgWriter.Write<double>(&image[0], width, inpReader.PhaseCentreRA(), inpReader.PhaseCentreDec(), inpReader.PixelSizeX(), inpReader.PixelSizeY());
	
	FitsWriter modelWriter(modelImageName);
	modelWriter.Write<double>(&model[0], width, inpReader.PhaseCentreRA(), inpReader.PhaseCentreDec(), inpReader.PixelSizeX(), inpReader.PixelSizeY());
}
