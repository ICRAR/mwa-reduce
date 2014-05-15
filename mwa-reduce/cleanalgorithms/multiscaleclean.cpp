#include "multiscaleclean.h"
#include "joinedclean.h"

#include "../lane.h"
#include "../fftresampler.h"
#include "../fftconvolver.h"

#include <boost/thread/thread.hpp>

#include <iostream>

template<typename ImageSetType>
void MultiScaleClean<ImageSetType>::ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain)
{
	if(this->_stopOnNegativeComponent)
		this->_allowNegativeComponents = true;
	_originalWidth = width;
	_originalHeight = height;
	
	_originalPsfs = &psfImages;
	_dataImageOriginal = &dataImage;
	_modelImage = &modelImage;
	
	_minScale = 4.0 * _beamSize / _pixelSizeX;
	double currentScale = _startScale;
	reachedStopGain = true;
	size_t iterationsAtStart = this->_iterationNumber;
	do {
		double nextScale = currentScale*0.5;
		std::cout << "Current scale: " << (currentScale*_pixelSizeX*180.0*60.0/M_PI) << " arcmin (" << round(currentScale) << " pixels).\n";
		bool reachedStopGainOnCurrentScale, canCleanFurther;
		size_t iterationsAtStartOfThisScale = this->_iterationNumber;
		executeMajorIterationForScale(currentScale, nextScale, reachedStopGainOnCurrentScale, canCleanFurther);
		std::cout << "Finished scale " << (currentScale*_pixelSizeX*180.0*60.0/M_PI) << " arcmin (" << round(currentScale) << " pixels), " << (this->_iterationNumber-iterationsAtStartOfThisScale) << " iterations performed, progressing to next scale.\n";
		if(canCleanFurther)
			reachedStopGain = reachedStopGain && reachedStopGainOnCurrentScale;
		currentScale = nextScale;
	} while(this->_iterationNumber < this->_maxIter && currentScale >= _minScale);
	
	if(this->_iterationNumber < this->_maxIter)
	{
		std::cout << "Minimum scale reached, finishing cleaning with delta functions.\n";
		// Finish off with normal clean
		JoinedClean<ImageSetType> joinedClean;
		joinedClean.CopyConfigFrom(*this);
		joinedClean.SetIterationNumber(this->IterationNumber());
		bool reachedStopGainOnCurrentScale;
		joinedClean.ExecuteMajorIteration(dataImage, modelImage, psfImages, width, height, reachedStopGainOnCurrentScale);
		reachedStopGain = reachedStopGainOnCurrentScale;
		this->_iterationNumber = joinedClean.IterationNumber();
	}
	if(reachedStopGain)
		std::cout << "Multiscale clean finished this major iteration after " << this->_iterationNumber-iterationsAtStart << " iterations; continuing after prediction/inversion round.\n";
	else
		std::cout << "Multiscale clean performed " << this->_iterationNumber-iterationsAtStart << " in this iteration. Absolute stopping criteria reached.\n";
}

template<typename ImageSetType>
void MultiScaleClean<ImageSetType>::executeMajorIterationForScale(double currentScale, double nextScale, bool& reachedStopGain, bool& canCleanFurther)
{
	ImageBufferAllocator<double>& allocator = *_dataImageOriginal->Allocator();
	// Scale down the image so that the "scale size" is 10 pixels
	// 10 pixels is so that the peak finding is still reasonably accurate
	_rescaledWidth = size_t(ceil(double(_originalWidth)*(10.0/currentScale))),
	_rescaledHeight = size_t(ceil(double(_originalHeight)*(10.0/currentScale)));
	double rescaleFactor = double(_rescaledWidth) / _originalWidth;
	if(_rescaledWidth > _originalWidth || _rescaledHeight > _originalHeight)
	{
		_rescaledWidth = _originalWidth;
		_rescaledHeight = _originalHeight;
		rescaleFactor = 1.0;
	}
	
	ImageSetType
		largeScaleImage(_rescaledWidth*_rescaledHeight, *_dataImageOriginal),
		nextScaleImage(_rescaledWidth*_rescaledHeight, *_dataImageOriginal),
		currentScaleModel(_rescaledWidth*_rescaledHeight, *_modelImage);
	std::vector<double*> scaledPsfs(_originalPsfs->size());
	_dataImageLargeScale = &largeScaleImage;
	_dataImageNextScale = &nextScaleImage;
	_scaledPsfs = &scaledPsfs;
	
	canCleanFurther = true;
	size_t cpuCount = this->_threadCount;
	
	double thresholdBias = pow(_thresholdBiasLevel, log2(currentScale));
	std::cout << "Threshold bias: " << thresholdBias << '\n';
	double oldSubtractionGain = this->_subtractionGain;
	this->_subtractionGain *= sqrt(thresholdBias);
	
	// Fill the large and next scale images with the rescaled images
	FFTResampler imageResampler(_originalWidth, _originalHeight, _rescaledWidth, _rescaledHeight, cpuCount, false);
	imageResampler.Start();
	for(size_t i=0; i!=_dataImageOriginal->ImageCount(); ++i)
		imageResampler.AddTask(_dataImageOriginal->GetImage(i), largeScaleImage.GetImage(i));
	for(size_t i=0; i!=_originalPsfs->size(); ++i)
	{
		scaledPsfs[i] = allocator.Allocate(_rescaledWidth*_rescaledHeight);
		imageResampler.AddTask((*_originalPsfs)[i], scaledPsfs[i]);
	}
	imageResampler.Finish();
	
	for(size_t i=0; i!=_dataImageOriginal->ImageCount(); ++i)
		memcpy(nextScaleImage.GetImage(i), largeScaleImage.GetImage(i), sizeof(double)*_rescaledWidth*_rescaledHeight);
	for(size_t i=0; i!=currentScaleModel.ImageCount(); ++i)
		memset(currentScaleModel.GetImage(i), 0, sizeof(double)*_rescaledWidth*_rescaledHeight);
	
	// Convolve the large and next scale images and the psfs
	double* kernelImage = allocator.Allocate(_rescaledWidth * _rescaledHeight);
	ao::uvector<double> shape;
	size_t kernelSize;
	MakeShapeFunction(currentScale * rescaleFactor, shape, kernelSize);
	memset(kernelImage, 0, sizeof(double) * _rescaledWidth * _rescaledHeight);
	FFTConvolver::PrepareKernel(kernelImage, _rescaledWidth, _rescaledHeight, shape.data(), kernelSize);
	for(size_t i=0; i!=largeScaleImage.ImageCount(); ++i)
		FFTConvolver::ConvolveSameSize(largeScaleImage.GetImage(i), kernelImage, _rescaledWidth, _rescaledHeight);
	for(size_t i=0; i!=scaledPsfs.size(); ++i) {
		FFTConvolver::ConvolveSameSize(scaledPsfs[i], kernelImage, _rescaledWidth, _rescaledHeight);
		FFTConvolver::ConvolveSameSize(scaledPsfs[i], kernelImage, _rescaledWidth, _rescaledHeight);
	}
	
	MakeShapeFunction(nextScale * rescaleFactor, shape, kernelSize);
	memset(kernelImage, 0, sizeof(double) * _rescaledWidth * _rescaledHeight);
	FFTConvolver::PrepareKernel(kernelImage, _rescaledWidth, _rescaledHeight, shape.data(), kernelSize);
	for(size_t i=0; i!=nextScaleImage.ImageCount(); ++i)
		FFTConvolver::ConvolveSameSize(nextScaleImage.GetImage(i), kernelImage, _rescaledWidth, _rescaledHeight);
	
	//FitsWriter writer;
	//writer.SetImageDimensions(_rescaledWidth, _rescaledHeight);
	//writer.Write("multiscale-kernel.fits", kernelImage);
	//writer.Write("multiscale-current-scale.fits", largeScaleImage.GetImage(0));
	//writer.Write("multiscale-next-scale.fits", nextScaleImage.GetImage(0));
	//writer.Write("multiscale-psf.fits", scaledPsfs[0]);
	//writer.SetImageDimensions(kernelSize, kernelSize);
	//writer.Write("multiscale-shape.fits", shape.data());
		
	size_t componentX=0, componentY=0;
	findPeak(componentX, componentY);
	std::cout << "Initial peak: " << peakDescription(largeScaleImage, componentX, componentY, rescaleFactor) << '\n';
	
	size_t peakIndex = componentX + componentY*_rescaledWidth;
	double peakNormalized = _dataImageLargeScale->JoinedValueNormalized(peakIndex) * rescaleFactor * rescaleFactor;
	double firstThreshold = this->_threshold, stopGainThreshold = peakNormalized*(1.0-this->_stopGain)/thresholdBias;
	std::cout << "Scale-adjusted threshold: " << firstThreshold*thresholdBias << ", major iteration stops at " << stopGainThreshold*thresholdBias << '\n';
	if(stopGainThreshold > firstThreshold)
	{
		firstThreshold = stopGainThreshold;
		std::cout << "Next major iteration for this scale at: " << stopGainThreshold << '\n';
	}
	else if(this->_stopGain != 1.0) {
		std::cout << "Major iteration threshold reached global threshold of " << this->_threshold << " for this scale.\n";
	}

	std::vector<ao::lane<CleanTask>*> taskLanes(cpuCount);
	std::vector<ao::lane<CleanResult>*> resultLanes(cpuCount);
	boost::thread_group threadGroup;
	for(size_t i=0; i!=cpuCount; ++i)
	{
		taskLanes[i] = new ao::lane<CleanTask>(1);
		resultLanes[i] = new ao::lane<CleanResult>(1);
		CleanThreadData cleanThreadData;
		cleanThreadData.startY = (_rescaledHeight*i)/cpuCount;
		cleanThreadData.endY = _rescaledHeight*(i+1)/cpuCount;
		cleanThreadData.parent = this;
		threadGroup.add_thread(new boost::thread(&MultiScaleClean<ImageSetType>::cleanThreadFunc, this, &*taskLanes[i], &*resultLanes[i], cleanThreadData));
	}
	
	while(fabs(peakNormalized) > firstThreshold*thresholdBias && this->_iterationNumber < this->_maxIter && !(largeScaleImage.IsComponentNegative(peakIndex) && this->_stopOnNegativeComponent))
	{
		if(this->_iterationNumber <= 10 ||
			(this->_iterationNumber <= 100 && this->_iterationNumber % 10 == 0) ||
			(this->_iterationNumber <= 1000 && this->_iterationNumber % 100 == 0) ||
			this->_iterationNumber % 1000 == 0)
			std::cout << "Iteration " << this->_iterationNumber << ": " << peakDescription(largeScaleImage, componentX, componentY, rescaleFactor) << '\n';
		
		CleanTask task;
		task.cleanCompX = componentX;
		task.cleanCompY = componentY;
		task.peak = largeScaleImage.Get(peakIndex);
		for(size_t i=0; i!=cpuCount; ++i)
			taskLanes[i]->write(task);
		
		currentScaleModel.AddComponent(largeScaleImage, peakIndex, this->_subtractionGain * rescaleFactor * rescaleFactor);
		
		double peakUnnormalized = 0.0;
		for(size_t i=0; i!=cpuCount; ++i)
		{
			CleanResult result;
			resultLanes[i]->read(result);
			if(result.peakLevelUnnormalized >= peakUnnormalized && std::isfinite(result.peakLevelUnnormalized))
			{
				peakUnnormalized = result.peakLevelUnnormalized;
				componentX = result.nextPeakX;
				componentY = result.nextPeakY;
			}
		}
		if(peakUnnormalized == 0.0)
		{
			std::cout << "No more components found at current scale: continuing to next scale.\n";
			canCleanFurther = false;
			break;
		}
		peakIndex = componentX + componentY*_rescaledWidth;
		peakNormalized = largeScaleImage.JoinedValueNormalized(peakIndex) * rescaleFactor * rescaleFactor;
		
		++this->_iterationNumber;
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
	reachedStopGain = fabs(peakNormalized) <= stopGainThreshold*thresholdBias;
	
	for(size_t i=0; i!=scaledPsfs.size(); ++i)
		allocator.Free(scaledPsfs[i]);
	
	MakeShapeFunction(currentScale * rescaleFactor, shape, kernelSize);
	memset(kernelImage, 0, sizeof(double) * _rescaledWidth * _rescaledHeight);
	FFTConvolver::PrepareKernel(kernelImage, _rescaledWidth, _rescaledHeight, shape.data(), kernelSize);
	
	double* convolvedModel = allocator.Allocate(_originalWidth * _originalHeight);
	double* preparedPsf = allocator.Allocate(_originalWidth * _originalHeight);
	FFTResampler modelResampler(_rescaledWidth, _rescaledHeight, _originalWidth, _originalHeight, cpuCount, false);
	for(size_t i=0; i!=currentScaleModel.ImageCount(); ++i)
	{
		FFTConvolver::ConvolveSameSize(currentScaleModel.GetImage(i), kernelImage, _rescaledWidth, _rescaledHeight);
		//writer.SetImageDimensions(_rescaledWidth, _rescaledHeight);
		//writer.Write("multiscale-model-small.fits", currentScaleModel.GetImage(i));
		modelResampler.RunSingle(currentScaleModel.GetImage(i), convolvedModel);
		double *modelPtr = convolvedModel;
		double *dest = _modelImage->GetImage(i), *destEnd = dest + _originalWidth*_originalHeight;
		while(dest != destEnd) {
			*dest += *modelPtr;
			++modelPtr;
			++dest;
		}
		FFTConvolver::PrepareForConvolution(preparedPsf, (*_originalPsfs)[_modelImage->PSFIndex(i)], _originalWidth, _originalHeight);
		FFTConvolver::ConvolveSameSize(convolvedModel, preparedPsf, _originalWidth, _originalHeight);
		dest = _dataImageOriginal->GetImage(i);
		destEnd = dest + _originalWidth*_originalHeight;
		modelPtr = convolvedModel;
		while(dest != destEnd) {
			*dest -= *modelPtr;
			++modelPtr;
			++dest;
		}
	}
	allocator.Free(preparedPsf);
	allocator.Free(convolvedModel);
	allocator.Free(kernelImage);
	
	this->_subtractionGain = oldSubtractionGain;
}

template<typename ImageSetType>
void MultiScaleClean<ImageSetType>::findPeak(size_t& x, size_t& y, size_t startY, size_t stopY) const
{
	double peakMax = std::numeric_limits<double>::min();
	const size_t lastIndex = _rescaledWidth*_rescaledHeight;
	size_t peakIndex = lastIndex;
	
	const double scaleBias = scaleBiasFunction(1.0, 2.0);
	const size_t
		horBorderSize = floor(_rescaledWidth*this->CleanBorderRatio()),
		verBorderSize = floor(_rescaledHeight*this->CleanBorderRatio());
	size_t xiStart = horBorderSize, xiEnd = _rescaledWidth - horBorderSize;
	size_t yiStart = std::max(startY, verBorderSize), yiEnd = std::min(stopY, _rescaledHeight - verBorderSize);
	if(xiEnd < xiStart) xiEnd = xiStart;
	if(yiEnd < yiStart) yiEnd = yiStart;
	for(size_t yi=yiStart; yi!=yiEnd; ++yi)
	{
		size_t index=yi*_rescaledWidth + xiStart;
		for(size_t xi=xiStart; xi!=xiEnd; ++xi)
		{
			double value = _dataImageLargeScale->JoinedValue(index);
			if(std::isfinite(value) && fabs(value) > peakMax)
			{
				double valueNextScale = _dataImageNextScale->JoinedValue(index);
				if(value < 0.0)
				{
					value = -value;
					valueNextScale = -valueNextScale;
				}
				if(std::isfinite(valueNextScale) && value > valueNextScale * scaleBias)
				{
					peakIndex = index;
					peakMax = value;
				}
			}
			++index;
		}
	}
	if(peakIndex == lastIndex)
	{
		x = _rescaledWidth;
		y = _rescaledHeight;
	}
	else {
		x = peakIndex % _rescaledWidth;
		y = peakIndex / _rescaledWidth;
	}
}

template<typename ImageSetType>
void MultiScaleClean<ImageSetType>::cleanThreadFunc(ao::lane<CleanTask> *taskLane, ao::lane<CleanResult> *resultLane, CleanThreadData cleanData)
{
	CleanTask task;
	// This initialization is not really necessary, but gcc warns about possible uninitialized values otherwise
	task.peak = ImageSetType::Value::Zero();
	ImageSetType& imageSet = *cleanData.parent->_dataImageLargeScale;
	ImageSetType& nextScaleSet = *cleanData.parent->_dataImageNextScale;
	const std::vector<double*>& psfs = *cleanData.parent->_scaledPsfs;
	while(taskLane->read(task))
	{
		for(size_t i=0; i!=imageSet.ImageCount(); ++i)
		{
			subtractImage(imageSet.GetImage(i), psfs[ImageSetType::PSFIndex(i)], task.cleanCompX, task.cleanCompY, this->_subtractionGain * task.peak.GetValue(i), cleanData.startY, cleanData.endY);
		}
		for(size_t i=0; i!=nextScaleSet.ImageCount(); ++i)
		{
			subtractImage(nextScaleSet.GetImage(i), psfs[ImageSetType::PSFIndex(i)], task.cleanCompX, task.cleanCompY, this->_subtractionGain * task.peak.GetValue(i), cleanData.startY, cleanData.endY);
		}
		
		CleanResult result;
		findPeak(result.nextPeakX, result.nextPeakY, cleanData.startY, cleanData.endY);
		if(result.nextPeakX == _rescaledWidth)
			result.peakLevelUnnormalized = std::numeric_limits<double>::quiet_NaN();
		else
			result.peakLevelUnnormalized = imageSet.AbsJoinedValue(result.nextPeakX + result.nextPeakY*_rescaledWidth);
		
		resultLane->write(result);
	}
}

template<typename ImageSetType>
std::string MultiScaleClean<ImageSetType>::peakDescription(const ImageSetType& image, size_t x, size_t y, double rescaleFactor)
{
	std::ostringstream str;
	size_t index = x + y*_rescaledWidth;
	double peak = image.JoinedValueNormalized(index) * rescaleFactor * rescaleFactor;
	str << peak << " Jy at " << x << "," << y;
	return str.str();
}

template class MultiScaleClean<clean_algorithms::PolarizedImageSet<2>>;
template class MultiScaleClean<clean_algorithms::PolarizedImageSet<4>>;
template class MultiScaleClean<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<2>>>;
template class MultiScaleClean<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<4>>>;
template class MultiScaleClean<clean_algorithms::SingleImageSet>;
