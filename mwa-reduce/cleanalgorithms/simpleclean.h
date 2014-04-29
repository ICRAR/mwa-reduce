#ifndef SIMPLE_CLEAN_H
#define SIMPLE_CLEAN_H

#include <string>
#include <cmath>
#include <limits>

#include "cleanalgorithm.h"
#include "imageset.h"

namespace ao {
	template<typename T> class lane;
}

class SimpleClean : public TypedCleanAlgorithm<clean_algorithms::SingleImageSet>
{
	public:
#ifdef __AVX__
		template<bool AllowNegativeComponent>
		static double FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, double borderRatio);
#endif
		
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, double borderRatio)
		{
			double peakMax = std::numeric_limits<double>::min();
			size_t peakIndex = width * height;
			
			const size_t horBorderSize = round(width*borderRatio), verBorderSize = round(height*borderRatio);
			const size_t xiStart = horBorderSize, xiEnd = width - horBorderSize;
			const size_t yiStart = verBorderSize, yiEnd = height - verBorderSize;
	
			for(size_t yi=yiStart; yi!=yiEnd; ++yi)
			{
				size_t index = yi*width + xiStart;;
				for(size_t xi=xiStart; xi!=xiEnd; ++xi)
				{
					double value = image[index];
					if(std::isfinite(value))
					{
						if(allowNegativeComponents) value = std::fabs(value);
						if(value > peakMax)
						{
							peakIndex = index;
							peakMax = std::fabs(value);
						}
					}
					++value;
				}
			}
			if(peakIndex == width * height)
			{
				x = width; y = height;
				return std::numeric_limits<double>::quiet_NaN();
			}
			else {
				x = peakIndex % width;
				y = peakIndex / width;
				return image[x + y*width];
			}
		}

		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, const class AreaSet &cleanAreas);

#ifdef __AVX__
		static double PartialFindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			double peakLevel;
			if(allowNegativeComponents)
				peakLevel = FindPeakAVX<true>(image + (width * startY), width, endY-startY, x, y, borderRatio);
			else
				peakLevel = FindPeakAVX<false>(image + (width * startY), width, endY-startY, x, y, borderRatio);
			y += startY;
			return peakLevel;
		}
		static double PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			return PartialFindPeakAVX(image, width, height, x, y, allowNegativeComponents, startY, endY, borderRatio);
		}
#else
		static double PartialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			return PartialFindPeakSimple(image, width, height, x, y, allowNegativeComponents, startY, endY, borderRatio);
		}
#endif

		static double PartialFindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			double peakLevel = FindPeak(image + (width * startY), width, endY-startY, x, y, allowNegativeComponents, borderRatio);
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
		
    virtual void ExecuteMajorIteration(ImageSet& dataImage, ImageSet& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain)
		{
			ExecuteMajorIteration(dataImage.GetImage(0), modelImage.GetImage(0), psfImages[0], width, height, reachedStopGain);
		}
		
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
