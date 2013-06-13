#ifndef CLEAN_ALGORITHM_H
#define CLEAN_ALGORITHM_H

#include <string>
#include <cmath>

class CleanAlgorithm
{
	public:
		CleanAlgorithm();
		
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents)
		{
			double peakMax = fabs(*image);
			const double *imgIter = image, *endPtr = image + width * height;
			size_t peakIndex = 0;
			size_t index = 0;
			while(!std::isfinite(peakMax) && imgIter!=endPtr)
			{
				peakMax = fabs(*imgIter);
				++imgIter;
				++index;
				++peakIndex;
			}
			for(const double *i=imgIter; i!=endPtr; ++i)
			{
				double value = *i;
				if(std::isfinite(value))
				{
					if(allowNegativeComponents) value = fabs(value);
					if(value > peakMax)
					{
						peakIndex = index;
						peakMax = fabs(*i);
					}
				}
				++index;
			}
			x = peakIndex % width;
			y = peakIndex / width;
			return image[x + y*width];
		}

		static void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
		
		void ExecuteMajorIteration(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height);
		
		void SetMaxNIter(size_t nIter) { _maxIter = nIter; }

	private:
		double _threshold, _subtractionGain;
		size_t _maxIter;
};

#endif
