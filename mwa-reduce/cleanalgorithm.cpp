#include "cleanalgorithm.h"
#include "imagecoordinates.h"
#include "modelsource.h"
#include "model.h"

#include <iostream>
#include <emmintrin.h>

CleanAlgorithm::CleanAlgorithm() :
	_threshold(0.0),
	_subtractionGain(0.1),
	_maxIter(500)
{
}

void CleanAlgorithm::SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor)
{
	size_t startX, startY, endX, endY;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
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
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = y + height/2;
	if(endY > height) endY = height;
	
	__m128d factor2 = _mm_set_pd(factor, factor);
	for(size_t ypos = startY; ypos != endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			// I've SSE-ified this, but it didn't improve speed at all :-/
			// (Compiler probably already did it)
			//*imageIter = *imageIter - (*psfIter * factor);
			//*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			_mm_storeu_pd(imageIter, _mm_sub_pd(_mm_loadu_pd(imageIter), _mm_mul_pd(_mm_loadu_pd(psfIter), factor2)));
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}


void CleanAlgorithm::ExecuteMajorIteration(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height)
{
	bool allowNegativeComponents = false;
	
	size_t componentX, componentY;
	double peak = CleanAlgorithm::FindPeak(dataImage, width, height, componentX, componentY, allowNegativeComponents);
	std::cout << "Initial peak: " << peak << '\n';
	size_t iterationNumber = 0;
	while(fabs(peak) > _threshold && iterationNumber < _maxIter)
	{
		if(iterationNumber % 10 == 0)
			std::cout << "Iteration " << iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		CleanAlgorithm::SubtractImage(dataImage, psfImage, width, height, componentX, componentY, _subtractionGain * peak);
		modelImage[componentX + componentY*width] += _subtractionGain * peak;
		
		peak = CleanAlgorithm::FindPeak(dataImage, width, height, componentX, componentY, allowNegativeComponents);
		++iterationNumber;
	}
	std::cout << "Stopped on peak " << peak << '\n';
}

void CleanAlgorithm::GetModelFromImage(Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double spectralIndex, double refFreq)
{
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			double value = image[y*width + x];
			if(value != 0.0 && std::isfinite(value))
			{
				long double l, m;
				ImageCoordinates::XYToLM<long double>(x, y, pixelSizeX, pixelSizeY, width, height, l, m);
			
				ModelSource source;
				long double ra, dec;
				ImageCoordinates::LMToRaDec<long double>(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
				std::stringstream nameStr;
				nameStr << "component" << model.SourceCount();
				source.SetName(nameStr.str());
				source.SetBrightness(SourceSDFWithSI<long double>(value, spectralIndex, refFreq));
				source.SetPosRA(ra);
				source.SetPosDec(dec);
				model.AddSource(source);
			}
		}
	}
}
