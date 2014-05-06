#include "fftconvolver.h"

#include "uvector.h"

#include <fftw3.h>

#include <complex>
#include <stdexcept>

void FFTConvolver::Convolve(double* image, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize)
{
	ao::uvector<double> scaledKernel(imgWidth * imgHeight, 0.0);
	PrepareKernel(scaledKernel.data(), imgWidth, imgHeight, kernel, kernelSize);
	ConvolveSameSize(image, scaledKernel.data(), imgWidth, imgHeight);
}

void FFTConvolver::PrepareKernel(double* dest, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize)
{
	if(kernelSize > imgWidth || kernelSize > imgHeight)
		throw std::runtime_error("Kernel size > image dimension");
	const double* kernelIter = kernel;
	for(size_t y=0; y!=kernelSize/2; ++y)
	{
		size_t destY = imgHeight - kernelSize/2 + y;
		size_t firstX = imgWidth - kernelSize/2;
		double *destIter = &dest[destY * imgWidth + firstX];
		for(size_t x=0; x!=kernelSize/2; ++x)
		{
			*destIter = *kernelIter;
			++kernelIter;
			++destIter;
		}
		destIter = &dest[destY * imgWidth];
		for(size_t x=kernelSize/2; x!=kernelSize; ++x)
		{
			*destIter = *kernelIter;
			++kernelIter;
			++destIter;
		}
	}
	for(size_t y=kernelSize/2; y!=kernelSize; ++y)
	{
		size_t firstX = imgWidth - kernelSize/2;
		double *destIter = &dest[firstX + (y-kernelSize/2)*imgWidth];
		for(size_t x=0; x!=kernelSize/2; ++x)
		{
			*destIter = *kernelIter;
			++kernelIter;
			++destIter;
		}
		destIter = &dest[(y-kernelSize/2)*imgWidth];
		for(size_t x=kernelSize/2; x!=kernelSize; ++x)
		{
			*destIter = *kernelIter;
			++kernelIter;
			++destIter;
		}
	}
}

void FFTConvolver::PrepareForConvolution(double* dest, const double* source, size_t imgWidth, size_t imgHeight)
{
	const double* sourceIter = source;
	for(size_t y=0; y!=imgHeight/2; ++y)
	{
		size_t destY = imgHeight - imgHeight/2 + y;
		size_t firstX = imgWidth - imgWidth/2;
		double *destIter = &dest[destY * imgWidth + firstX];
		for(size_t x=0; x!=imgWidth/2; ++x)
		{
			*destIter = *sourceIter;
			++sourceIter;
			++destIter;
		}
		destIter = &dest[destY * imgWidth];
		for(size_t x=imgWidth/2; x!=imgWidth; ++x)
		{
			*destIter = *sourceIter;
			++sourceIter;
			++destIter;
		}
	}
	for(size_t y=imgHeight/2; y!=imgHeight; ++y)
	{
		size_t firstX = imgWidth - imgWidth/2;
		double *destIter = &dest[firstX + (y-imgHeight/2)*imgWidth];
		for(size_t x=0; x!=imgWidth/2; ++x)
		{
			*destIter = *sourceIter;
			++sourceIter;
			++destIter;
		}
		destIter = &dest[(y-imgHeight/2)*imgWidth];
		for(size_t x=imgWidth/2; x!=imgWidth; ++x)
		{
			*destIter = *sourceIter;
			++sourceIter;
			++destIter;
		}
	}
}

void FFTConvolver::ConvolveSameSize(double* image, const double* kernel, size_t imgWidth, size_t imgHeight)
{
	const size_t imgSize = imgWidth * imgHeight;
	double* tempData = reinterpret_cast<double*>(fftw_malloc(imgSize * sizeof(double)));
	fftw_complex* fftImageData = reinterpret_cast<fftw_complex*>(fftw_malloc(imgSize * sizeof(fftw_complex)));
	fftw_complex* fftKernelData = reinterpret_cast<fftw_complex*>(fftw_malloc(imgSize * sizeof(fftw_complex)));
	fftw_plan inToFPlan = fftw_plan_dft_r2c_2d(imgHeight, imgWidth, tempData, fftImageData, FFTW_ESTIMATE);
	fftw_plan fToOutPlan = fftw_plan_dft_c2r_2d(imgHeight, imgWidth, fftImageData, tempData, FFTW_ESTIMATE);
	
	memcpy(tempData, image, imgSize * sizeof(double));
	fftw_execute_dft_r2c(inToFPlan, tempData, reinterpret_cast<fftw_complex*>(fftImageData));
	
	memcpy(tempData, kernel, imgSize * sizeof(double));
	fftw_execute_dft_r2c(inToFPlan, tempData, reinterpret_cast<fftw_complex*>(fftKernelData));
	
	double fact = 1.0/imgSize;
	for(size_t i=0; i!=imgSize; ++i)
		reinterpret_cast<std::complex<double>*>(fftImageData)[i] *= fact * reinterpret_cast<std::complex<double>*>(fftKernelData)[i];
		
	fftw_execute_dft_c2r(fToOutPlan, reinterpret_cast<fftw_complex*>(fftImageData), tempData);
	memcpy(image, tempData, imgSize * sizeof(double));
		
	fftw_free(fftImageData);
	fftw_free(fftKernelData);
	fftw_free(tempData);
	
	fftw_destroy_plan(inToFPlan);
	fftw_destroy_plan(fToOutPlan);
}

