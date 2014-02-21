#ifndef JOINED_POL_CLEAN_H
#define JOINED_POL_CLEAN_H

#include "cleanalgorithm.h"
#include "simpleclean.h"

#include "../imagebufferallocator.h"
#include "../cachedimageset.h"

namespace ao {
	template<typename T> class lane;
}

class JoinedPolClean : public CleanAlgorithm
{
public:
	class ImageSet {
	public:
		double *xx, *xyr, *xyi, *yy;
		
		ImageSet(size_t size, ImageBufferAllocator<double>& allocator) :
			xx(allocator.Allocate(size)),
			xyr(allocator.Allocate(size)),
			xyi(allocator.Allocate(size)),
			yy(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		~ImageSet()
		{
			_allocator->Free(xx);
			_allocator->Free(xyr);
			_allocator->Free(xyi);
			_allocator->Free(yy);
		}
		
		void Load(CachedImageSet& set)
		{
			set.Load(xx, PolarizationEnum::XX, false);
			set.Load(xyr, PolarizationEnum::XY, false);
			set.Load(xyi, PolarizationEnum::XY, true);
			set.Load(yy, PolarizationEnum::YY, false);
		}
		
		void Store(CachedImageSet& set) const
		{
			set.Store(xx, PolarizationEnum::XX, false);
			set.Store(xyr, PolarizationEnum::XY, false);
			set.Store(xyi, PolarizationEnum::XY, true);
			set.Store(yy, PolarizationEnum::YY, false);
		}
		
		double* Get(size_t index)
		{
			double* vals[4] = { xx, xyr, xyi, yy };
			return vals[index];
		}
		
		double SquaredSum(size_t index) const
		{
			return
				xx[index]*xx[index] +
				xyr[index]*xyr[index] + xyi[index]*xyi[index] +
				yy[index]*yy[index];
		}
		
		bool IsComponentNegative(size_t index) const
		{
			return xx[index]<0.0 || yy[index]<0.0;
		}
		
		void AddComponent(const ImageSet& source, size_t index, double factor)
		{
			xx[index] += source.xx[index] * factor;
			xyr[index] += source.xyr[index] * factor;
			xyi[index] += source.xyi[index] * factor;
			yy[index] += source.yy[index] * factor;
		}
		
	private:
		ImageBufferAllocator<double> *_allocator;
	};
		
	void ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain);
private:
	size_t _width, _height;
	
	struct CleanTask
	{
		size_t cleanCompX, cleanCompY;
		double peakXX, peakXYr, peakXYi, peakYY;
	};
	struct CleanResult
	{
		CleanResult() : nextPeakX(0), nextPeakY(0), peakLevelSquared(0.0)
		{ }
		size_t nextPeakX, nextPeakY;
		double peakLevelSquared;
	};
	struct CleanThreadData
	{
		size_t startY, endY;
		ImageSet* dataImage;
		const double* psfImage;
	};

	void findPeak(const ImageSet& image, size_t& x, size_t& y) const
	{
		findPeak(image, x, y, 0, _height);
	}
	void findPeak(const ImageSet& image, size_t& x, size_t& y, size_t startY, size_t stopY) const;
	
	std::string peakDescription(const ImageSet& image, size_t& x, size_t& y);
	void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _width, _height, psf, _width, _height, x, y, factor, startY, endY);
	}
};

#endif
