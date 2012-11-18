
#include "btpimager.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <stdint.h>

#define FFTW_PLAN_COUNT 64

#ifdef __SSE__
	#include <xmmintrin.h>
	#include <emmintrin.h>
	#include <smmintrin.h>
	#define DO_VECTORIZE true
#else
	#define DO_VECTORIZE false
#endif

#define USE_DFT false
#define TWO_STAGE_FFT true

using std::size_t;

BTPImager::BTPImager(size_t imageCount, size_t imgSize, NumType pixelScale) :
	_fftInputs(0),
	_fftOutputs(0),
	_fftwPlansC2R(0),
	_fftwPlansR2C(0),
	_fftwPlanSizesC2R(0),
	_fftwPlanSizesR2C(0),
	_imageData(new ImageNum*[imageCount]),
	_lookupSqrtLMTerm(0),
	_isInitialized(false),
	_pixelScale(pixelScale),
	_imgSize(imgSize),
	_imageCount(imageCount),
	_weightCounters(new NumType[_imageCount]),
	_overallMaxUVDist(0.0)
{
	if(_imgSize%4 != 0)
	{
		_imgSize = _imgSize + (3 - ((_imgSize-1)%4));
		std::cout << "WARNING: Image size was not an exact multiply of 4. Will scale image size up to " << _imgSize << '\n';
	}
	
	for(size_t image=0;image!=_imageCount;++image)
	{
		_imageData[image] = (ImageNum*) fftw_malloc(sizeof(ImageNum) * _imgSize * _imgSize);
		_weightCounters[image] = 0.0;
		
		for(size_t i=0;i!=_imgSize*_imgSize;++i)
		{
			_imageData[image][i] = 0.0;
		}
	}
	
	_lookupSqrtLMTerm = new FastNum[_imgSize*_imgSize];
	size_t lookupI = 0;
	NumType midX = (NumType) _imgSize / 2.0, midY = (NumType) _imgSize / 2.0;
	NumType scaleLM = _pixelScale;
	for(size_t y=0;y!=_imgSize;++y)
	{
		for(size_t x=0;x!=_imgSize;++x)
		{
			NumType l = ((NumType) x - midX) * scaleLM;
			NumType m = (midY - (NumType) y) * scaleLM;
			NumType r = m*m + l*l;
			if(r < 1.0)
				_lookupSqrtLMTerm[lookupI] = sqrt(1.0-r)-1.0;
			else
				_lookupSqrtLMTerm[lookupI] = -1.0;
			++lookupI;
		}
	}
}

BTPImager::~BTPImager()
{
	deinitialize();
	
	for(size_t image=0;image!=_imageCount;++image)
	{
		fftw_free(_imageData[image]);
	}
	
	delete[] _imageData;
	delete[] _weightCounters;
	delete[] _lookupSqrtLMTerm;
}

void BTPImager::Initialize(NumType smallestUVDistTimesLambda, NumType largestUVDistTimesLambda, NumType highestFrequency, NumType frequencyStep, size_t channelCount)
{
	deinitialize();
	
	_highestFrequency = highestFrequency;
	_frequencyStep = frequencyStep;
	_channelCount = channelCount;
	_sampleDist = 2*((size_t) round(highestFrequency/frequencyStep));

	if(fmod(highestFrequency, frequencyStep) > frequencyStep*0.0001 && !USE_DFT)
	{
		std::cout << "WARNING: For now, to be able to perform exact gridding without DFT, your channels frequencies need to be\n"
			"an exact multiply of the frequency step! (offset: " << fmod(highestFrequency, frequencyStep) << " Hz)\n";
	}
	_startChannel = (size_t) round(highestFrequency/_frequencyStep - (channelCount-1));
	_minLambda = frequencyToWavelength(_highestFrequency);
	NumType uvDist = smallestUVDistTimesLambda / _minLambda;
	_largestFFTSize = std::max((size_t) (_sampleDist/(_pixelScale * (2.0*uvDist))), 2*_sampleDist);
	_largestFFTSize = _largestFFTSize + (15 - ((_largestFFTSize-1)%16));
	
	_fftInputs = new FftwComplex*[_imageCount];
	_fftOutputs = new FftwNum*[_imageCount];
	for(size_t image=0;image!=_imageCount;++image)
	{
		_fftInputs[image] = (FftwComplex*) fftw_malloc(sizeof(FftwComplex) * (_largestFFTSize/2+1));
		_fftOutputs[image] = (FftwNum*) fftw_malloc(sizeof(FftwNum) * _largestFFTSize);
	}
	
	//fftw_import_wisdom_from_filename("fftw_wisdom.txt");
	size_t currentFFTSize = _largestFFTSize;
	_fftwPlansC2R = new fftw_plan[FFTW_PLAN_COUNT];
	_fftwPlanSizesC2R = new size_t[FFTW_PLAN_COUNT];
	_fftwPlansR2C = new fftw_plan[FFTW_PLAN_COUNT];
	_fftwPlanSizesR2C = new size_t[FFTW_PLAN_COUNT];
	
	for(size_t i=0;i<FFTW_PLAN_COUNT;++i)
	{
		_fftwPlansC2R[i] = fftw_plan_dft_c2r_1d(currentFFTSize, _fftInputs[0], _fftOutputs[0], FFTW_ESTIMATE); // FFTW_MEASURE ?
		_fftwPlanSizesC2R[i] = currentFFTSize;
		_fftwPlansR2C[i] = fftw_plan_dft_r2c_1d(currentFFTSize, _fftOutputs[0], _fftInputs[0], FFTW_ESTIMATE); // FFTW_MEASURE ?
		_fftwPlanSizesR2C[i] = currentFFTSize;
		currentFFTSize = size_t(currentFFTSize/sqrt(sqrt(2)));
		currentFFTSize = currentFFTSize + (15 - ((currentFFTSize-1)%16));
	}
	
	//fftw_export_wisdom_to_filename("fftw_wisdom.txt");
	std::cout << "Largest FFT: " << _largestFFTSize << '\n';
	
	NumType maxUVDist = largestUVDistTimesLambda / frequencyToWavelength(highestFrequency - 0.5*frequencyStep*channelCount);
	if(maxUVDist > _overallMaxUVDist)
		_overallMaxUVDist = maxUVDist;
	
	_isInitialized = true;
}

void BTPImager::deinitialize()
{
	if(_isInitialized)
	{
		for(size_t image=0;image!=_imageCount;++image)
		{
			fftw_free(_fftInputs[image]);
			fftw_free(_fftOutputs[image]);
		}
		for(size_t i=0;i!=FFTW_PLAN_COUNT;++i)
		{
			fftw_destroy_plan(_fftwPlansC2R[i]);
			fftw_destroy_plan(_fftwPlansR2C[i]);
		}
		delete[] _fftwPlansC2R;
		delete[] _fftwPlanSizesC2R;
		delete[] _fftwPlansR2C;
		delete[] _fftwPlanSizesR2C;
		
		delete[] _fftInputs;
		delete[] _fftOutputs;
		
		_isInitialized = false;
	}
}

void BTPImager::AddTimestep(size_t imageIndex, NumType uTimesLambda, NumType vTimesLambda, NumType wTimesLambda, const std::complex<float> *data, NumType zenithDistance, NumType paralacticAngle, NumType weight)
{
	//if(imageIndex != 0 || _weightCounters[0]!=0.0) return;
	NumType phi = atan2(vTimesLambda, uTimesLambda);
	NumType uvDist = sqrt(uTimesLambda*uTimesLambda + vTimesLambda*vTimesLambda) / _minLambda;
	
	// Steps to be taken:
	// 1. Create a 1D array with the data in it (in 'u' dir) and otherwise zerod.
	//    This represents a line through the uvw plane. The highest non-zero samples
	//    have uv-distance 2*|uvw|. Channels are not regridded. Therefore, the distance in
	//    samples is 2 * (lowestFrequency / frequencyStep + (channelCount-1))
	//    (Another two times larger than necessary to prevent border issues).
	//    
	// 2. Fourier transform this (FFT 1D)
	
	// 3. Stretch, rotate and make it fill the plane
	//    - Stretch with 1/|uvw| and wide field corrections
	//    - Rotate with phi
	// 4. Add to output
	
	size_t minimumFFTSize = 2*_sampleDist;
	//size_t minimumFFTSize = std::max((size_t) (_sampleDist/(_pixelScale * (2.0*uvDist))), 2*_sampleDist);
	//if(minimumFFTSize / _sampleDist > 1.0)
	//	minimumFFTSize = _sampleDist * 1.0;
	
	size_t fftSize;
	fftw_plan fftwPlan;
	if(USE_DFT)
		fftSize = minimumFFTSize;
	else
		getFFTWPlanComplexToReal(minimumFFTSize, fftSize, fftwPlan);
	
	// Calculate the transformation terms
	NumType cosPhi = cos(phi), sinPhi = sin(phi);
	NumType mid = (NumType) _imgSize / 2.0;
	NumType transformGen = (2.0*uvDist / _sampleDist) * fftSize;
	NumType transformX = cosPhi * transformGen;
	// Negative rotation (thus positive sin sign)
	NumType transformY = sinPhi * transformGen;
	ImageNum *destPtr = _imageData[imageIndex];
	NumType scaleLM = _pixelScale;
	FastNum *lookupSqrtLMTermPtr = _lookupSqrtLMTerm;
	// Wide field correction factors
	NumType tanZcosChi = tan(zenithDistance) * cos(paralacticAngle);
	NumType tanZsinChi = tan(zenithDistance) * sin(paralacticAngle);
	
	double *fftOut = _fftOutputs[imageIndex];
	size_t fftCentre = fftSize / 2;
	
	// Calculate extends (maximum and minimum index into FFT array)
	NumType
		maxM = mid * scaleLM,
		minM = (mid - (NumType) (_imgSize-1)) * scaleLM,
		maxMW = std::max(maxM, maxM + tanZcosChi),
		minMW = std::min(minM, minM + tanZcosChi),
		maxL = ((NumType) (_imgSize-1) - mid) * scaleLM,
		minL = -mid * scaleLM,
		maxLW = std::max(maxL, maxL - tanZsinChi),
		minLW = std::min(minL, minL - tanZsinChi);
	double
		maxSrcXDb = round((maxLW + maxMW) * transformGen) + fftCentre + 1,
		minSrcXDb = round((minLW + minMW) * transformGen) + fftCentre;
	fftw_complex *fftInp = _fftInputs[imageIndex];
			
	//std::cout << "Extends: " << minSrcXDb << " - " << maxSrcXDb << " of " << fftSize << '\n';
	if(USE_DFT)
	{
		size_t maxSrcX, minSrcX;
		if(minSrcXDb <= 0.0 || maxSrcXDb >= fftSize)
		{
			minSrcX = 0;
			maxSrcX = fftSize;
		} else {
			minSrcX = size_t(minSrcXDb);
			maxSrcX = size_t(maxSrcXDb);
		}
		
		for(size_t x = minSrcX; x != maxSrcX; ++x)
			fftOut[x] = 0.0;
			
		NumType lowestFrequency = _highestFrequency - (_frequencyStep * (_channelCount-1));
		double fftFactor = sqrt(_channelCount*2) / fftSize;
		for(size_t ch = 0; ch != _channelCount; ++ch)
		{
			double
				valueReal = data[ch].real() * fftFactor,
				valueImag = data[ch].imag() * fftFactor,
				expChFactor = 2.0 * M_PI * (lowestFrequency / _frequencyStep + ch) / fftSize;
			for(size_t x = minSrcX; x != maxSrcX; ++x)
			{
				double angle = expChFactor * (int) (x - fftCentre);
				double sinAngle, cosAngle;
				sincos(angle, &sinAngle, &cosAngle);
				fftOut[x] += valueReal * cosAngle - valueImag * sinAngle;
			}
		}
	} else {
		// USE FFT
		
		for(size_t i=0;i!=(fftSize/2+1);++i) {
			fftInp[i][0] = 0.0;
			fftInp[i][1] = 0.0;
		}
		
		// fftw gives unnormalized results; have to divide by sqrt n.
		// Have to divide by two, because we add the series twice (non-conjugate and conjugate)
		// but will do that when normalizing images.
		NumType fftFactor = 1.0; //1.0 / sqrt(_channelCount);
		NumType mulFactor = fftFactor;// * _channelCount*sqrt(2) / (_sampleDist/(_pixelScale * (2.0*uvDist)));
		const size_t startChannel = _startChannel, channelCount = _channelCount;
		for(size_t ch=0;ch!=channelCount;++ch)
		{
			fftInp[(startChannel + ch)][0] = data->real() * mulFactor;
			fftInp[(startChannel + ch)][1] = data->imag() * mulFactor;
			++data;
		}
		
		fftw_execute_dft_c2r(fftwPlan, fftInp, fftOut);
		
		for(size_t ch=0;ch!=fftCentre;++ch)
			std::swap(fftOut[ch+fftCentre], fftOut[ch]);
	}
	
	double superSamplingFactor = 4.0;
	double minResampleFactor = superSamplingFactor * _sampleDist/(_pixelScale * (2.0*uvDist) * fftSize);
	if(TWO_STAGE_FFT && minResampleFactor > 1.0)
	{
		// Our FFT output array contains fewer samples per distance than the
		// image plane. To avoid interpolation errors, the part of the array that is used is
		// scaled up by sinc convolution resampling.
		size_t leftX, rightX;
		const double secStageOversize = 2.0;
		if(minSrcXDb <= 0.0 || (maxSrcXDb-fftCentre)*secStageOversize+fftCentre >= fftSize)
		{
			leftX = 0;
			rightX = fftSize;
			std::cout << "Full size 2nd stage FFT!\n";
		} else {
			double thisCentre = round((maxSrcXDb + minSrcXDb) * 0.5);
			leftX = size_t(round((minSrcXDb-thisCentre)*secStageOversize+thisCentre));
			rightX = size_t(round((maxSrcXDb-thisCentre)*secStageOversize+thisCentre+1));
		}
		
		size_t resampleR2CSize;
		fftw_plan resamplePlanR2C;
		getFFTWPlanRealToComplex(rightX-leftX, resampleR2CSize, resamplePlanR2C);
		size_t extraSpace = resampleR2CSize - rightX + leftX;
		size_t extraLeftSpace = extraSpace / 2;
		if(leftX < extraLeftSpace)
		{
			extraLeftSpace = leftX;
			leftX = 0;
			std::cout << "Extension did not fit on left side!\n";
		} else {
			leftX -= extraLeftSpace;
		}
		if(rightX + (extraSpace - extraLeftSpace) >= fftSize) {
			rightX = fftSize;
			std::cout << "Extension did not fit!\n";
		} else
			rightX += (extraSpace - extraLeftSpace);
		
		for(size_t x=leftX;x!=rightX;++x)
			fftOut[x-leftX] = fftOut[x];
		for(size_t x=rightX-leftX;x!=resampleR2CSize;++x)
			fftOut[x] = 0.0;
		double minResampleDestSize = resampleR2CSize*minResampleFactor;
		size_t resampleC2RSize;
		fftw_plan resamplePlanC2R;
		getFFTWPlanComplexToReal(minResampleDestSize, resampleC2RSize, resamplePlanC2R);
		for(size_t x=resampleR2CSize/2+1;x!=resampleC2RSize/2+1;++x)
		{
			fftInp[x][0] = 0.0;
			fftInp[x][1] = 0.0;
		}
		fftw_execute_dft_r2c(resamplePlanR2C, fftOut, fftInp);
		fftw_execute_dft_c2r(resamplePlanC2R, fftInp, fftOut);
		
		double
			fftNormFactorA = 1.0 / sqrt(resampleR2CSize),
			fftNormFactorB = 1.0 / sqrt((double) resampleR2CSize),
			fftNormFactor = fftNormFactorA * fftNormFactorB;
		for(size_t x=0;x!=resampleC2RSize;++x)
			fftOut[x] *= fftNormFactor;
		
		// Now change the transform parameters to accomodate the new size of the fft array.
		double resampleFactor = (double) resampleC2RSize / (double) resampleR2CSize;
		size_t tempCentre = (size_t) round(((double) fftCentre - (double) leftX) * resampleFactor);
		//std::cout << "Scaled " << leftX << '-' << rightX << " to " << resampleC2RSize << ", " << resampleFactor << "x center:" << fftCentre << "->" << tempCentre << "fftSize=" << (resampleC2RSize + tempCentre) << '\n';
		fftCentre = tempCentre;
		fftSize = resampleC2RSize + tempCentre;
		transformX *= resampleFactor;
		transformY *= resampleFactor;
	} else {
		std::cout << "Did not supersample.\n";
	}
	
	if(!DO_VECTORIZE) {
		size_t srcXLLimit = (size_t) -1, srcXULimit = 0;
		for(size_t y=0;y!=_imgSize;++y)
		{
			NumType m = (mid - (NumType) y) * scaleLM;
			
			for(size_t x=0;x!=_imgSize;++x)
			{
				NumType sqrtlmTerm = *lookupSqrtLMTermPtr;
				NumType mw = m - sqrtlmTerm * tanZcosChi;
				NumType yrTransformed = mw * transformY;
			
				NumType l = ((NumType) x - mid) * scaleLM;
				NumType lw = l + sqrtlmTerm * tanZsinChi;
				NumType srcX = lw * transformX + yrTransformed;
				
				if(false) { //DEBUG
					size_t srcXIndexLimit = (size_t) round(srcX) + fftCentre;
					if(srcXIndexLimit > srcXULimit) srcXULimit = srcXIndexLimit;
					if(srcXIndexLimit < srcXLLimit) srcXLLimit = srcXIndexLimit;
				}
				
				size_t srcXIndex = ((size_t) round(srcX) + fftCentre) % fftSize;
				*destPtr += fftOut[srcXIndex];
				++destPtr;
				++lookupSqrtLMTermPtr;
			}
		}
		//std::cout << "LLimit: " << srcXLLimit << " ULimit: " << srcXULimit << '\n';
	} else {
		/** VECTORIZED VERSION */
		const __m128 tanZcosChi_ps = _mm_set1_ps(tanZcosChi);
		const __m128 tanZsinChi_ps = _mm_set1_ps(tanZsinChi);
		const __m128 transformX_ps = _mm_set1_ps(transformX);
		const __m128 transformY_ps = _mm_set1_ps(transformY);
		const __m128 fftCentre_ps = _mm_set1_ps(fftCentre);
		const __m128 fftSize_ps = _mm_set1_ps(fftSize);
		const __m128 scaleLM_ps = _mm_set1_ps(scaleLM);
		const __m128 mid_ps = _mm_set1_ps(mid);
		for(size_t y=0;y!=_imgSize;++y)
		{
			const __m128 m = _mm_set1_ps((mid - (FastNum) y) * scaleLM);
			
			for(size_t x=0;x!=_imgSize;x+=4)
			{
				__m128 sqrtlmTerm = _mm_load_ps(lookupSqrtLMTermPtr);
				__m128 mw = _mm_sub_ps(m, _mm_mul_ps(sqrtlmTerm, tanZcosChi_ps));
				__m128 yrTransformed = _mm_mul_ps(mw, transformY_ps);
			
				__m128 l = _mm_mul_ps(_mm_sub_ps(_mm_set_ps((FastNum) (x+3), (FastNum) (x+2), (FastNum) (x+1), (FastNum) x), mid_ps), scaleLM_ps);
				__m128 lw = _mm_add_ps(l, _mm_mul_ps(sqrtlmTerm, tanZsinChi_ps));
				__m128 srcX = _mm_add_ps(_mm_mul_ps(lw, transformX_ps), yrTransformed);
				srcX = _mm_add_ps(srcX, fftCentre_ps);
				// Calculate floor(num/fftSize)*fftSize (i.e. num - (num % fftSize))
				__m128 modTerm = _mm_mul_ps(_mm_floor_ps(_mm_div_ps(srcX, fftSize_ps)), fftSize_ps);
				// Calculate round(srcX - modTerm)
				__m128i srcXi = _mm_cvtps_epi32(_mm_sub_ps(srcX, modTerm));
				union { __m128i m; uint32_t i[4]; } srcXpi;
				_mm_store_si128(&srcXpi.m, srcXi);
				
				for(size_t count=0;count!=4;++count)
				{
					uint32_t srcXIndex = srcXpi.i[count];
					*destPtr += fftOut[srcXIndex];
					++destPtr;
				}
				
				lookupSqrtLMTermPtr += 4;
			}
		}
	}
	// Times two because conjugate and non-conjugate were added.
	_weightCounters[imageIndex] += weight * _channelCount * 2.0;
}

void BTPImager::GetIntermediateResult(ImageNum *imageData)
{
	size_t totalSize = _imgSize*_imgSize;
	memcpy(imageData, _imageData[0], sizeof(ImageNum)*totalSize);
	NumType weight = _weightCounters[0];
	for(size_t i=1;i!=_imageCount;++i)
	{
		ImageNum *data = _imageData[i];
		for(size_t yx=0;yx!=totalSize;++yx)
		{
			imageData[yx] += data[yx];
		}
		weight += _weightCounters[i];
	}
	NumType factor = (weight!=0.0) ? 1.0 / weight : 1.0;
	NumType mid = (NumType) _imgSize / 2.0;
	
	size_t yx=0;
	for(size_t y=0;y!=_imgSize;++y)
	{
		for(size_t x=0;x!=_imgSize;++x)
		{
			NumType l = ((NumType) x - mid) * _pixelScale;
			NumType m = (mid - (NumType) y) * _pixelScale;
			imageData[yx] = factor * imageData[yx] * sqrt(1.0 - l*l - m*m);
			++yx;
		}
	}
}

void BTPImager::getFFTWPlanComplexToReal(size_t minimumFFTSize, size_t &actualFFTSize, fftw_plan &fftwPlan)
{
	actualFFTSize = _fftwPlanSizesC2R[0];
	fftwPlan = _fftwPlansC2R[0];
	for(size_t i=1;i!=FFTW_PLAN_COUNT;++i)
	{
		if(_fftwPlanSizesC2R[i] >= minimumFFTSize && _fftwPlanSizesC2R[i] < actualFFTSize)
		{
			actualFFTSize = _fftwPlanSizesC2R[i];
			fftwPlan = _fftwPlansC2R[i];
		}
	}
	if(actualFFTSize < minimumFFTSize)
		throw std::runtime_error("Something went wrong; could not find appropriate fft plan...");
}

void BTPImager::getFFTWPlanRealToComplex(size_t minimumFFTSize, size_t &actualFFTSize, fftw_plan &fftwPlan)
{
	actualFFTSize = _fftwPlanSizesR2C[0];
	fftwPlan = _fftwPlansR2C[0];
	for(size_t i=1;i!=FFTW_PLAN_COUNT;++i)
	{
		if(_fftwPlanSizesR2C[i] >= minimumFFTSize && _fftwPlanSizesR2C[i] < actualFFTSize)
		{
			actualFFTSize = _fftwPlanSizesR2C[i];
			fftwPlan = _fftwPlansR2C[i];
		}
	}
	if(actualFFTSize < minimumFFTSize)
		throw std::runtime_error("Something went wrong; could not find appropriate fft plan...");
}
