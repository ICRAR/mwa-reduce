#include "joinedclean.h"

#include "../lane.h"

#include <boost/thread/thread.hpp>
#include <emmintrin.h>
#include <immintrin.h>

template<typename ImageSetType>
void JoinedClean<ImageSetType>::ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain)
{
	if(_stopOnNegativeComponent)
		_allowNegativeComponents = true;
	_width = width;
	_height = height;
	
	size_t componentX=0, componentY=0;
	findPeak(dataImage, componentX, componentY);
	std::cout << "Initial peak: " << peakDescription(dataImage, componentX, componentY) << '\n';
	
	size_t peakIndex = componentX + componentY*_width;
	double peakNormalized = dataImage.JoinedValueNormalized(peakIndex);
	double firstThreshold = _threshold, stopGainThreshold = peakNormalized*(1.0-_stopGain);
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
		cleanThreadData.dataImage = &dataImage;
		cleanThreadData.psfImages = psfImages;
		cleanThreadData.startY = (height*i)/cpuCount;
		cleanThreadData.endY = height*(i+1)/cpuCount;
		threadGroup.add_thread(new boost::thread(&JoinedClean::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	
	while(peakNormalized > firstThreshold && _iterationNumber < _maxIter && !(dataImage.IsComponentNegative(peakIndex) && _stopOnNegativeComponent))
	{
		if(_iterationNumber <= 10 ||
			(_iterationNumber <= 100 && _iterationNumber % 10 == 0) ||
			(_iterationNumber <= 1000 && _iterationNumber % 100 == 0) ||
			_iterationNumber % 1000 == 0)
			std::cout << "Iteration " << _iterationNumber << ": " << peakDescription(dataImage, componentX, componentY) << '\n';
		
		CleanTask task;
		task.cleanCompX = componentX;
		task.cleanCompY = componentY;
		task.peak = dataImage.Get(peakIndex);
		for(size_t i=0; i!=cpuCount; ++i)
			taskLanes[i]->write(task);
		
		modelImage.AddComponent(dataImage, peakIndex, _subtractionGain);
		
		double peakUnnormalized = 0.0;
		for(size_t i=0; i!=cpuCount; ++i)
		{
			CleanResult result;
			resultLanes[i]->read(result);
			if(result.peakLevelUnnormalized >= peakUnnormalized)
			{
				peakUnnormalized = result.peakLevelUnnormalized;
				componentX = result.nextPeakX;
				componentY = result.nextPeakY;
			}
		}
		peakIndex = componentX + componentY*_width;
		peakNormalized = dataImage.JoinedValueNormalized(peakIndex);
		
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
	std::cout << "Stopped on peak " << peakNormalized << '\n';
	reachedStopGain = peakNormalized < stopGainThreshold;
}

template<typename ImageSetType>
void JoinedClean<ImageSetType>::findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY) const
{
	double peakMax = std::numeric_limits<double>::min();
	size_t peakIndex = 0;
	const size_t lastIndex = _width*_height;
	
	for(size_t index=0; index!=lastIndex; ++index)
	{
		double value = image.JoinedValue(index);
		if(std::isfinite(value))
		{
			if(value > peakMax)
			{
				peakIndex = index;
				peakMax = value;
			}
		}
	}
	x = peakIndex % _width;
	y = peakIndex / _width;
}

template<typename ImageSetType>
void JoinedClean<ImageSetType>::cleanThreadFunc(ao::lane<CleanTask> *taskLane, ao::lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	// This initialization is not really necessary, but gcc warns about possible uninitialized values otherwise
	task.peak = ImageSetType::Value::Zero();
	while(taskLane->read(task))
	{
		for(size_t i=0; i!=cleanData.dataImage->ImageCount(); ++i)
		{
			subtractImage(cleanData.dataImage->GetImage(i), cleanData.psfImages[ImageSetType::PSFIndex(i)], task.cleanCompX, task.cleanCompY, _subtractionGain * task.peak.GetValue(i), cleanData.startY, cleanData.endY);
		}
		
		CleanResult result;
		findPeak(*cleanData.dataImage, result.nextPeakX, result.nextPeakY, cleanData.startY, cleanData.endY);
		result.peakLevelUnnormalized = cleanData.dataImage->JoinedValue(result.nextPeakX + result.nextPeakY*_width);
		
		resultLane->write(result);
	}
}

template<typename ImageSetType>
std::string JoinedClean<ImageSetType>::peakDescription(const ImageSetType& image, size_t& x, size_t& y)
{
	std::ostringstream str;
	size_t index = x + y*_width;
	double peak = image.JoinedValueNormalized(index);
	str << peak << " Jy at " << x << "," << y;
	return str.str();
}

template class JoinedClean<joined_pol_clean::SingleImageSet>;
template class JoinedClean<joined_pol_clean::MultiImageSet>;