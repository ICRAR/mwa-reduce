#include "cleanalgorithm.h"
#include "imagecoordinates.h"
#include "modelsource.h"
#include "model.h"
#include "lane.h"
#include "areaset.h"

#include <boost/thread/thread.hpp>
#include <emmintrin.h>

#include <iostream>
#include <limits>

CleanAlgorithm::CleanAlgorithm() :
	_threshold(0.0),
	_subtractionGain(0.1),
	_stopGain(1.0),
	_maxIter(500),
	_iterationNumber(0),
	_allowNegativeComponents(false),
	_cleanAreas(0)
{
}

double CleanAlgorithm::partialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const class AreaSet &cleanAreas)
{
	double peakMax = std::numeric_limits<double>::min();
	size_t index = 0;
	const double *imgIter = &image[startY*width];
	x = 0; y = startY;
	for(size_t yi=startY; yi!=endY; ++yi)
	{
		for(size_t xi=0; xi!=width; ++xi)
		{
			double value = *imgIter;
			if(std::isfinite(value))
			{
				if(allowNegativeComponents) value = std::fabs(value);
				if(value > peakMax && cleanAreas.AllowCleaningInImage(xi, yi))
				{
					x = xi;
					y = yi;
					peakMax = std::fabs(*imgIter);
				}
			}
			++index;
			++imgIter;
		}
	}
	return image[x + y*width];
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

void CleanAlgorithm::partialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - width/2, offsetY = (int) y - height/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = x + width/2;
	if(endX > width) endX = width;
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + height/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * width + startX;
		const double *psfIter = psf + (ypos - offsetY) * width + startX - offsetX;
		for(size_t xpos = startX; xpos != endX; xpos+=2)
		{
			*imageIter = *imageIter - (*psfIter * factor);
			*(imageIter+1) = *(imageIter+1) - (*(psfIter+1) * factor);
			imageIter+=2;
			psfIter+=2;
		}
		if(!isAligned)
			*imageIter -= *psfIter * factor;
	}
}

void CleanAlgorithm::ExecuteMajorIterationST(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height)
{
	size_t componentX, componentY;
	double peak = FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents);
	std::cout << "Initial peak: " << peak << '\n';
	while(fabs(peak) > _threshold && _iterationNumber < _maxIter)
	{
		if(_iterationNumber % 10 == 0)
			std::cout << "Iteration " << _iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		SubtractImage(dataImage, psfImage, width, height, componentX, componentY, _subtractionGain * peak);
		modelImage[componentX + componentY*width] += _subtractionGain * peak;
		
		peak = FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents);
		++_iterationNumber;
	}
	std::cout << "Stopped on peak " << peak << '\n';
}

void CleanAlgorithm::ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain)
{
	size_t componentX, componentY;
	double peak = FindPeak(dataImage, width, height, componentX, componentY, _allowNegativeComponents);
	std::cout << "Initial peak: " << peak << '\n';
	double firstThreshold = _threshold, stopGainThreshold = fabs(peak*(1.0-_stopGain));
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		std::cout << "Next major iteration at: " << stopGainThreshold << '\n';
	}
	else if(_stopGain != 1.0) {
		std::cout << "Major iteration threshold reached global threshold of " << _threshold << ": final major iteration.\n";
	}

	size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	std::vector<lane<CleanTask>*> taskLanes(cpuCount);
	std::vector<lane<CleanResult>*> resultLanes(cpuCount);
	boost::thread_group threadGroup;
	for(size_t i=0; i!=cpuCount; ++i)
	{
		taskLanes[i] = new lane<CleanTask>(1);
		resultLanes[i] = new lane<CleanResult>(1);
		CleanThreadData cleanThreadData;
		cleanThreadData.width = width;
		cleanThreadData.height = height;
		cleanThreadData.dataImage = dataImage;
		cleanThreadData.psfImage = psfImage;
		cleanThreadData.startY = (height*i)/cpuCount;
		cleanThreadData.endY = height*(i+1)/cpuCount;
		threadGroup.add_thread(new boost::thread(&CleanAlgorithm::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	while(fabs(peak) > firstThreshold && _iterationNumber < _maxIter)
	{
		if(_iterationNumber % 10 == 0)
			std::cout << "Iteration " << _iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
		
		CleanTask task;
		task.cleanCompX = componentX;
		task.cleanCompY = componentY;
		task.peakLevel = peak;
		for(size_t i=0; i!=cpuCount; ++i)
			taskLanes[i]->write(task);
		
		modelImage[componentX + componentY*width] += _subtractionGain * peak;
		
		peak = 0.0;
		for(size_t i=0; i!=cpuCount; ++i)
		{
			CleanResult result;
			resultLanes[i]->read(result);
			if(fabs(result.peakLevel) >= fabs(peak))
			{
				peak = result.peakLevel;
				componentX = result.nextPeakX;
				componentY = result.nextPeakY;
			}
		}
		
		++_iterationNumber;
	}
	for(size_t i=0; i!=cpuCount; ++i)
		taskLanes[i]->write_end();
	threadGroup.join_all();
	for(size_t i=0; i!=cpuCount; ++i)
	{
		delete taskLanes[i];
		delete resultLanes[i];
	}
	std::cout << "Stopped on peak " << peak << '\n';
	reachedStopGain = fabs(peak) < stopGainThreshold;
}

void CleanAlgorithm::cleanThreadFunc(lane<CleanTask> *taskLane, lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	while(taskLane->read(task))
	{
		partialSubtractImage(cleanData.dataImage, cleanData.psfImage, cleanData.width, cleanData.height, task.cleanCompX, task.cleanCompY, _subtractionGain * task.peakLevel, cleanData.startY, cleanData.endY);
		
		CleanResult result;
		if(_cleanAreas == 0)
			result.peakLevel = partialFindPeak(cleanData.dataImage, cleanData.width, cleanData.height, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY);
		else
			result.peakLevel = partialFindPeak(cleanData.dataImage, cleanData.width, cleanData.height, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY, *_cleanAreas);
		
		resultLane->write(result);
	}
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
			
				ModelComponent component;
				long double ra, dec;
				ImageCoordinates::LMToRaDec<long double>(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
				std::stringstream nameStr;
				nameStr << "component" << model.SourceCount();
				component.SetSED(SpectralEnergyDistribution(value, refFreq, spectralIndex));
				component.SetPosRA(ra);
				component.SetPosDec(dec);
				
				ModelSource source;
				source.SetName(nameStr.str());
				source.AddComponent(component);
				model.AddSource(source);
			}
		}
	}
}

void CleanAlgorithm::PreparePSF(double* psf, size_t width, size_t height)
{
	double* endPtr = psf + width*height;
	while(psf != endPtr)
	{
		if(!std::isfinite(*psf)) *psf = 0.0;
		++psf;
	}
}

