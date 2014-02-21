#ifndef CLEAN_ALGORITHM_H
#define CLEAN_ALGORITHM_H

#include <string>
#include <cmath>

#include "../polarizationenum.h"

namespace ao {
	template<typename T> class lane;
}

class CleanAlgorithm
{
	public:
		virtual ~CleanAlgorithm() { }
		
		void SetMaxNIter(size_t nIter) { _maxIter = nIter; }
		
		void SetThreshold(double threshold) { _threshold = threshold; }
		
		void SetSubtractionGain(double gain) { _subtractionGain = gain; }
		
		void SetStopGain(double stopGain) { _stopGain = stopGain; }
		
		void SetAllowNegativeComponents(bool allowNegativeComponents) { _allowNegativeComponents = allowNegativeComponents; }
		
		void SetStopOnNegativeComponents(bool stopOnNegative) { _stopOnNegativeComponent = stopOnNegative; }
		
		void SetResizePSF(bool resizePSF) { _resizePSF = resizePSF; }
		
		size_t MaxNIter() const { return _maxIter; }
		double Threshold() const { return _threshold; }
		double SubtractionGain() const { return _subtractionGain; }
		double StopGain() const { return _stopGain; }
		bool AllowNegativeComponents() const { return _allowNegativeComponents; }
		bool StopOnNegativeComponents() const { return _allowNegativeComponents; }
		bool ResizePSF() const { return _resizePSF; }
		
		void SetCleanAreas(const class AreaSet& cleanAreas) { _cleanAreas = &cleanAreas; }
		
		size_t IterationNumber() const { return _iterationNumber; }
		
		
		static void ResizeImage(double* dest, size_t newWidth, size_t newHeight, const double* source, size_t width, size_t height);
		
		static void GetModelFromImage(class Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double spectralIndex, double refFreq, 
																	PolarizationEnum polarization = Polarization::StokesI);

		static void RemoveNaNsInPSF(double* psf, size_t width, size_t height);
		
		static void CalculateFastCleanPSFSize(size_t& psfWidth, size_t& psfHeight, size_t imageWidth, size_t imageHeight);
	protected:
		CleanAlgorithm();
		
		double _threshold, _subtractionGain, _stopGain;
		size_t _maxIter, _iterationNumber;
		bool _allowNegativeComponents, _stopOnNegativeComponent, _resizePSF;
		const class AreaSet *_cleanAreas;
};

#endif
