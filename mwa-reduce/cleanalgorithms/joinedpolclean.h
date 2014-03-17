#ifndef JOINED_POL_CLEAN_H
#define JOINED_POL_CLEAN_H

#include "cleanalgorithm.h"
#include "simpleclean.h"

#include "../imagebufferallocator.h"
#include "../cachedimageset.h"

namespace ao {
	template<typename T> class lane;
}

namespace joined_pol_clean {
		
	class SingleImageSet {
	public:
		double *xx, *xyr, *xyi, *yy;
		
		SingleImageSet(size_t size, ImageBufferAllocator<double>& allocator) :
			xx(allocator.Allocate(size)),
			xyr(allocator.Allocate(size)),
			xyi(allocator.Allocate(size)),
			yy(allocator.Allocate(size)),
			_allocator(&allocator)
		{
		}
		
		~SingleImageSet()
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
		
		double JoinedValue(size_t index) const
		{
			return SquaredSum(index);
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
		
		void AddComponent(const SingleImageSet& source, size_t index, double factor)
		{
			xx[index] += source.xx[index] * factor;
			xyr[index] += source.xyr[index] * factor;
			xyi[index] += source.xyi[index] * factor;
			yy[index] += source.yy[index] * factor;
		}
		
	private:
		ImageBufferAllocator<double> *_allocator;
	};
	
	class MultiImageSet {
	public:
		MultiImageSet(size_t imageSize, size_t count, ImageBufferAllocator<double>& allocator)
		{
			for(size_t i=0; i!=count; ++i)
			{
				_sets.push_back(new SingleImageSet(imageSize, allocator));
			}
		}
		
		~MultiImageSet()
		{
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				delete *i;
			}
		}
		
		void Load(CachedImageSet& set, size_t i)
		{
			_sets[i]->Load(set);
		}
		
		void Store(CachedImageSet& set, size_t i) const
		{
			_sets[i]->Store(set);
		}
		
		double JoinedValue(size_t index) const
		{
			double val = 0.0;
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				val += (*i)->SquaredSum(index);
			}
			return val;
		}
		
		bool IsComponentNegative(size_t index) const
		{
			for(std::vector<SingleImageSet*>::const_iterator i=_sets.begin(); i!=_sets.end(); ++i)
			{
				if((*i)->IsComponentNegative(index)) return true;
			}
			return false;
		}
		
		void AddComponent(const MultiImageSet& source, size_t index, double factor)
		{
			for(size_t i=0; i!=_sets.size(); ++i)
			{
				_sets[i]->AddComponent(*source._sets[i], index, factor);
			}
		}
		
	private:
		std::vector<SingleImageSet*> _sets;
	};	
}

template<typename ImageSetType = joined_pol_clean::SingleImageSet>
class JoinedPolClean : public CleanAlgorithm
{
public:
	void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain);
	typedef ImageSetType ImageSet;
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
		ImageSetType* dataImage;
		const double* psfImage;
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
