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
	_allowNegativeComponents(true),
	_stopOnNegativeComponent(false),
	_cleanAreas(0)
{
}

double CleanAlgorithm::PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const class AreaSet &cleanAreas)
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

void CleanAlgorithm::PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY)
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

void CleanAlgorithm::PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY)
{
	size_t startX, endX;
	int offsetX = (int) x - psfWidth/2, offsetY = (int) y - psfHeight/2;
	
	if(offsetX > 0)
		startX = offsetX;
	else
		startX = 0;
	
	if(offsetY > (int) startY)
		startY = offsetY;
	
	endX = std::min(x + psfWidth/2, imgWidth);
	
	bool isAligned = ((endX - startX) % 2) == 0;
	if(!isAligned) --endX;
	
	endY = std::min(y + psfHeight/2, endY);
	
	for(size_t ypos = startY; ypos < endY; ++ypos)
	{
		double *imageIter = image + ypos * imgWidth + startX;
		const double *psfIter = psf + (ypos - offsetY) * psfWidth + startX - offsetX;
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
	std::vector<double> resizedPsf;
	size_t psfWidth, psfHeight;
	if(_resizePSF)
	{
		CalculateFastCleanPSFSize(psfWidth, psfHeight, width, height);
		if(psfWidth != width || psfHeight != height)
		{
			resizedPsf.resize(psfWidth * psfHeight);
			ResizeImage(&resizedPsf[0], psfWidth, psfHeight, psfImage, width, height);
			psfImage = &resizedPsf[0];
		}
	}
	else {
		psfWidth = width;
		psfHeight = height;
	}
	if(_stopOnNegativeComponent)
		_allowNegativeComponents = true;
	
	size_t componentX=0, componentY=0;
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
	std::vector<ao::lane<CleanTask>*> taskLanes(cpuCount);
	std::vector<ao::lane<CleanResult>*> resultLanes(cpuCount);
	boost::thread_group threadGroup;
	for(size_t i=0; i!=cpuCount; ++i)
	{
		taskLanes[i] = new ao::lane<CleanTask>(1);
		resultLanes[i] = new ao::lane<CleanResult>(1);
		CleanThreadData cleanThreadData;
		cleanThreadData.imgWidth = width;
		cleanThreadData.imgHeight = height;
		cleanThreadData.dataImage = dataImage;
		cleanThreadData.psfWidth = psfWidth;
		cleanThreadData.psfHeight = psfHeight;
		cleanThreadData.psfImage = psfImage;
		cleanThreadData.startY = (height*i)/cpuCount;
		cleanThreadData.endY = height*(i+1)/cpuCount;
		threadGroup.add_thread(new boost::thread(&CleanAlgorithm::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	while(fabs(peak) > firstThreshold && _iterationNumber < _maxIter && (peak >= 0.0 || !_stopOnNegativeComponent))
	{
		if(
			(_iterationNumber <= 100 && _iterationNumber % 10 == 0) ||
			(_iterationNumber <= 1000 && _iterationNumber % 100 == 0) ||
			_iterationNumber % 1000 == 0)
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

void CleanAlgorithm::cleanThreadFunc(ao::lane<CleanTask> *taskLane, ao::lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	while(taskLane->read(task))
	{
		PartialSubtractImage(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, cleanData.psfImage, cleanData.psfWidth, cleanData.psfHeight, task.cleanCompX, task.cleanCompY, _subtractionGain * task.peakLevel, cleanData.startY, cleanData.endY);
		
		CleanResult result;
		if(_cleanAreas == 0)
			result.peakLevel = PartialFindPeak(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY);
		else
			result.peakLevel = PartialFindPeak(cleanData.dataImage, cleanData.imgWidth, cleanData.imgHeight, result.nextPeakX, result.nextPeakY, _allowNegativeComponents, cleanData.startY, cleanData.endY, *_cleanAreas);
		
		resultLane->write(result);
	}
}

void CleanAlgorithm::GetModelFromImage(Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double spectralIndex, double refFreq, PolarizationEnum polarization)
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
				component.SetSED(SpectralEnergyDistribution(value, refFreq, spectralIndex, polarization));
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

void CleanAlgorithm::ResizeImage(double* dest, size_t newWidth, size_t newHeight, const double* source, size_t width, size_t height)
{
	size_t srcStartX = (width - newWidth) / 2, srcStartY = (height - newHeight) / 2;
	for(size_t y=0; y!=newHeight; ++y)
	{
		double* destPtr = dest + y * newWidth;
		const double* srcPtr = source + (y + srcStartY) * width + srcStartX;
		memcpy(destPtr, srcPtr, newWidth * sizeof(double));
	}
}

void CleanAlgorithm::RemoveNaNsInPSF(double* psf, size_t width, size_t height)
{
	double* endPtr = psf + width*height;
	while(psf != endPtr)
	{
		if(!std::isfinite(*psf)) *psf = 0.0;
		++psf;
	}
}

void CleanAlgorithm::CalculateFastCleanPSFSize(size_t& psfWidth, size_t& psfHeight, size_t imageWidth, size_t imageHeight)
{
	// With 2048 x 2048, the subtraction is already so quick that it is not really required to make the psf smaller
	if(imageWidth <= 2048)
		psfWidth = imageWidth;
	else if(imageWidth <= 4096)
		psfWidth = 2048;
	else
		psfWidth = imageWidth / 2;
	
	if(imageHeight <= 2048)
		psfHeight = imageHeight;
	else if(imageHeight <= 4096)
		psfHeight = 2048;
	else
		psfHeight = imageHeight / 2;
}
