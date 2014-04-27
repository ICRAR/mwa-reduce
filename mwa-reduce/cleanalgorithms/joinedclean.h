#ifndef JOINED_CLEAN_H
#define JOINED_CLEAN_H

#include "cleanalgorithm.h"
#include "imageset.h"
#include "simpleclean.h"

namespace ao {
	template<typename T> class lane;
}

template<typename ImageSetType>
class JoinedClean : public TypedCleanAlgorithm<ImageSetType>
{
public:
	virtual void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain);
	
private:
	size_t _width, _height;
	
	struct CleanTask
	{
		size_t cleanCompX, cleanCompY;
		typename ImageSetType::Value peak;
	};
	struct CleanResult
	{
		CleanResult() : nextPeakX(0), nextPeakY(0), peakLevelUnnormalized(0.0)
		{ }
		size_t nextPeakX, nextPeakY;
		double peakLevelUnnormalized;
	};
	struct CleanThreadData
	{
		size_t startY, endY;
		ImageSetType* dataImage;
		std::vector<double*> psfImages;
	};

	void findPeak(const ImageSetType& image, size_t& x, size_t& y) const
	{
		findPeak(image, x, y, 0, _height);
	}
	void findPeak(const ImageSetType& image, size_t& x, size_t& y, size_t startY, size_t stopY) const;
	
	std::string peakDescription(const ImageSetType& image, size_t& x, size_t& y);
	void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
};

#endif
