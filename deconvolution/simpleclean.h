#ifndef SIMPLE_CLEAN_H
#define SIMPLE_CLEAN_H

#include <string>
#include <cmath>
#include <limits>

#include <boost/optional/optional.hpp>

#include "deconvolutionalgorithm.h"

#ifdef __SSE__
#define USE_INTRINSICS
#endif

namespace ao {
	template<typename T> class lane;
}

class SimpleClean
{
	public:
		SimpleClean() = delete;
		
		static boost::optional<double> FindPeakSimple(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
		
#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
		template<bool AllowNegativeComponent>
		static boost::optional<double> FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder);
		
		static boost::optional<double> FindPeakAVX(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder)
		{
			if(allowNegativeComponents)
				return FindPeakAVX<true>(image, width, height, x, y, startY, endY, horizontalBorder, verticalBorder);
			else
				return FindPeakAVX<false>(image, width, height, x, y, startY, endY, horizontalBorder, verticalBorder);
		}
#endif

		static boost::optional<double> FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, double borderRatio)
		{
			return FindPeak(image, width, height, x, y, allowNegativeComponents, startY, endY, round(width*borderRatio), round(height*borderRatio));
		}

		static boost::optional<double> FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, size_t horizontalBorder, size_t verticalBorder)
		{
#if defined __AVX__ && defined USE_INTRINSICS && !defined FORCE_NON_AVX
			return FindPeakAVX(image, width, height, x, y, allowNegativeComponents, startY, endY, horizontalBorder, verticalBorder);
#else
			return FindPeakSimple(image, width, height, x, y, allowNegativeComponents, startY, endY, horizontalBorder, verticalBorder);
#endif
		}

		static boost::optional<double> FindPeakWithMask(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, const bool* cleanMask);

		static boost::optional<double> FindPeakWithMask(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const bool* cleanMask, double borderRatio)
		{
			return FindPeakWithMask(image, width, height, x, y, allowNegativeComponents, startY, endY, cleanMask, round(width*borderRatio), round(height*borderRatio));
		}
		
		static boost::optional<double> FindPeakWithMask(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const bool* cleanMask, size_t horizontalBorder, size_t verticalBorder);
		
		static void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
		
		static void PartialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
		static void PartialSubtractImage(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
		
#if defined __AVX__ && defined USE_INTRINSICS
		static void PartialSubtractImageAVX(double *image, size_t imgWidth, size_t imgHeight, const double *psf, size_t psfWidth, size_t psfHeight, size_t x, size_t y, double factor, size_t startY, size_t endY);
#endif
		
	private:
};

#endif
