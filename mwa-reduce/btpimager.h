#ifndef BTPIMAGER_H
#define BTPIMAGER_H

#include <cstring>
#include <complex>

#include <fftw3.h>

class BTPImager
{
	public:
		typedef long double NumType;
		typedef double ImageNum;
		
		BTPImager(size_t imageCount, size_t imgSize, NumType pixelScale);
		~BTPImager();
		
		void Initialize(NumType smallestUVDistTimesLambda, NumType largestUVDistTimesLambda, NumType highestFrequency, NumType frequencyStep, std::size_t channelCount);
		
		void AddTimestep(size_t imageIndex, NumType uTimesLambda, NumType vTimesLambda, NumType wTimesLambda, const std::complex<float> *data, NumType zenithDistance, NumType paralacticAngle, NumType weight);
		
		void GetIntermediateResult(ImageNum *imageData);
		
		std::size_t ImageCount() const { return _imageCount; }
		std::size_t ImageSize() const { return _imgSize; }
		NumType PixelScale() const { return _pixelScale; }
		template<typename T>
		static T frequencyToWavelength(const T frequency)
		{
			return speedOfLight() / frequency; 
		}
		NumType OverallMaxUVDist() const { return _overallMaxUVDist; }
	private:
		typedef double FftwNum;
		typedef fftw_complex FftwComplex;
		typedef float FastNum;
		
		static long double speedOfLight()
		{
			return 299792458.0L;
		}
		void getFFTWPlanComplexToReal(std::size_t minimumFFTSize, std::size_t &actualFFTSize, fftw_plan &fftwPlan);
		void getFFTWPlanRealToComplex(std::size_t minimumFFTSize, std::size_t &actualFFTSize, fftw_plan &fftwPlan);
		void deinitialize();
		
		FftwComplex **_fftInputs;
		FftwNum **_fftOutputs;
		fftw_plan *_fftwPlansC2R, *_fftwPlansR2C;
		std::size_t *_fftwPlanSizesC2R, *_fftwPlanSizesR2C;
		ImageNum **_imageData;
		FastNum *_lookupSqrtLMTerm;
		bool _isInitialized;

		NumType _highestFrequency, _frequencyStep, _pixelScale, _minLambda;
		std::size_t _channelCount, _imgSize, _sampleDist, _largestFFTSize, _startChannel, _imageCount;
		NumType *_weightCounters;
		
		NumType _overallMaxUVDist;
};

#endif
