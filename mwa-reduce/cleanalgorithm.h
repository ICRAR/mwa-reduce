#ifndef CLEAN_ALGORITHM_H
#define CLEAN_ALGORITHM_H

#include <string>
#include <cmath>

template<typename T> class lane;

class CleanAlgorithm
{
	public:
		CleanAlgorithm();
		
		static double FindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents)
		{
			double peakMax = std::fabs(*image);
			const double *imgIter = image, *endPtr = image + width * height;
			size_t peakIndex = 0;
			size_t index = 0;
			while(!std::isfinite(peakMax) && imgIter!=endPtr)
			{
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

		static void SubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor);
		
		/**
		 * Single threaded implementation -- just for reference.
		 */
		void ExecuteMajorIterationST(double *dataImage, double *modelImage, const double *psfImage, size_t width, size_t height);
		
		void ExecuteMajorIteration(double* dataImage, double* modelImage, const double* psfImage, size_t width, size_t height, bool& reachedStopGain);
		
		void SetMaxNIter(size_t nIter) { _maxIter = nIter; }
		
		void SetThreshold(double threshold) { _threshold = threshold; }
		
		void SetSubtractionGain(double gain) { _subtractionGain = gain; }
		
		void SetStopGain(double stopGain) { _stopGain = stopGain; }
		
		void SetAllowNegativeComponents(bool allowNegativeComponents) { _allowNegativeComponents = allowNegativeComponents; }
		
		static void GetModelFromImage(class Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double spectralIndex, double refFreq);

		void SetCleanAreas(const class AreaSet& cleanAreas) { _cleanAreas = &cleanAreas; }
	private:
		struct CleanTask
		{
			size_t cleanCompX, cleanCompY;
			double peakLevel;
		};
		struct CleanResult
		{
			size_t nextPeakX, nextPeakY;
			double peakLevel;
		};
		struct CleanThreadData
		{
			size_t startY, endY;
			size_t width, height;
			double *dataImage;
			const double *psfImage;
		};
		void cleanThreadFunc(lane<CleanTask> *taskLane, lane<CleanResult> *resultLane, CleanThreadData cleanData);
		
		static void partialSubtractImage(double *image, const double *psf, size_t width, size_t height, size_t x, size_t y, double factor, size_t startY, size_t endY);
		static double partialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY)
		{
			double peakLevel = FindPeak(image + (width * startY), width, endY-startY, x, y, allowNegativeComponents);
			y += startY;
			return peakLevel;
		}
		static double partialFindPeak(const double *image, size_t width, size_t height, size_t &x, size_t &y, bool allowNegativeComponents, size_t startY, size_t endY, const class AreaSet &cleanAreas);
		
		double _threshold, _subtractionGain, _stopGain;
		size_t _maxIter, _iterationNumber;
		bool _allowNegativeComponents;
		const class AreaSet *_cleanAreas;
};

#endif
