#ifndef SIMPLE_CLEAN_H
#define SIMPLE_CLEAN_H

#include <string>
#include <cmath>

#include "cleanalgorithm.h"

namespace ao {
	template<typename T> class lane;
}

class SimpleClean : public CleanAlgorithm
{
	public:
#ifdef __AVX__
		template<bool AllowNegativeComponent>
		static double FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y);
#endif
		
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents)
		{
			double peakMax = std::fabs(*image);
			const double* imgIter = image;
			const double* const endPtr = image + width * height;
			size_t peakIndex = 0;
			size_t index = 0;
			
			while(!std::isfinite(peakMax) && imgIter!=endPtr)
			{
				if(allowNegativeComponents || *imgIter >= 0.0)
					peakMax = std::fabs(*imgIter);
				++imgIter;
				++index;
				++peakIndex;
			}
			for(const double *i=imgIter; i!=endPtr; ++i)
			{
				double value = *i;
				if(std::isfinite(value))
				{
					if(allowNegativeComponents) value = std::fabs(value);
					if(value > peakMax)
					{
						peakIndex = index;
						peakMax = std::fabs(*i);
					}
				}
				++index;
			}
			x = peakIndex % width;
			y = peakIndex / width;
			return image[x + y*width];
		}

		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, const class AreaSet &cleanAreas);

#ifdef __AVX__
		static double PartialFindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY)
		{
			double peakLevel;
			if(allowNegativeComponents)
				peakLevel = FindPeakAVX<true>(image + (width * startY), width, endY-startY, x, y);
			else
				peakLevel = FindPeakAVX<false>(image + (width * startY), width, endY-startY, x, y);
			y += startY;
			return peakLevel;
		}
		static double PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY)
		{
			return PartialFindPeakAVX(image, width, height, x, y, allowNegativeComponents, startY, endY);
		}
#else
		static double PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY)
		{
			return PartialFindPeakSimple(image, width, height, x, y, allowNegativeComponents, startY, endY);
		}
#endif

		static double PartialFindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY)
		{
			double peakLevel = FindPeak(image + (width * startY), width, endY-startY, x, y, allowNegativeComponents);
			y += startY;
			return peakLevel;
		}
		
		static double PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const class AreaSet &cleanAreas);
		
		static void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
		
		static void PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
		static void PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
#ifdef __AVX__
		static void PartialSubtractImageAVX(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
#endif
		
		/**
		 * Single threaded implementation -- just for reference.
		 */
		void ExecuteMajorIterationST(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height);
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain);
	private:
		struct CleanTask
		{
			size_t cleanCompX, cleanCompY;
			double peakLevel;
		};
		struct CleanResult
		{
			CleanResult() : nextPeakX(0), nextPeakY(0), peakLevel(0.0)
			{ }
			size_t nextPeakX, nextPeakY;
			double peakLevel;
		};
		struct CleanThreadData
		{
			size_t startY, endY;
			double *dataImage;
			size_t imgWidth, imgHeight;
			const double *psfImage;
			size_t psfWidth, psfHeight;
		};
		void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
};

#endif
