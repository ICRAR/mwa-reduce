#include "fftconvolver.h"

#include "fftwmanager.h"
#include "uvector.h"

#include <fftw3.h>

#include <complex>
#include <stdexcept>

void FFTConvolver::Convolve(FFTWManager& fftw, double* image, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize)
{
	ao::uvector<double> scaledKernel(imgWidth * imgHeight, 0.0);
	PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel, kernelSize);
	ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight);
}

void FFTConvolver::ReverseAndConvolve(class FFTWManager& fftw, double* image, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize)
{
	ao::uvector<double> scaledKernel(imgWidth * imgHeight, 0.0);
	
	PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel, kernelSize);
	ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight);
}

void FFTConvolver::PrepareSmallKernel(double* dest, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize)
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

void FFTConvolver::PrepareKernel(double* dest, const double* source, size_t imgWidth, size_t imgHeight)
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

void FFTConvolver::ConvolveSameSize(FFTWManager& fftw, double* image, const double* kernel, size_t imgWidth, size_t imgHeight)
{
	const size_t imgSize = imgWidth * imgHeight;
	const size_t complexSize = (imgWidth/2+1) * imgHeight;
	double* tempData = reinterpret_cast<double*>(fftw_malloc(imgSize * sizeof(double)));
	fftw_complex* fftImageData = reinterpret_cast<fftw_complex*>(fftw_malloc(complexSize * sizeof(fftw_complex)));
	fftw_complex* fftKernelData = reinterpret_cast<fftw_complex*>(fftw_malloc(complexSize * sizeof(fftw_complex)));
	
	std::unique_lock<std::mutex> lock(fftw.Mutex());
	fftw_plan inToFPlan = fftw_plan_dft_r2c_2d(imgHeight, imgWidth, tempData, fftImageData, FFTW_ESTIMATE);
	fftw_plan fToOutPlan = fftw_plan_dft_c2r_2d(imgHeight, imgWidth, fftImageData, tempData, FFTW_ESTIMATE);
	lock.unlock();
	
	memcpy(tempData, image, imgSize * sizeof(double));
	fftw_execute_dft_r2c(inToFPlan, tempData, fftImageData);
	
	memcpy(tempData, kernel, imgSize * sizeof(double));
	fftw_execute_dft_r2c(inToFPlan, tempData, fftKernelData);
	
	double fact = 1.0/imgSize;
	for(size_t i=0; i!=complexSize; ++i)
		reinterpret_cast<std::complex<double>*>(fftImageData)[i] *= fact * reinterpret_cast<std::complex<double>*>(fftKernelData)[i];
		
	fftw_execute_dft_c2r(fToOutPlan, reinterpret_cast<fftw_complex*>(fftImageData), tempData);
	memcpy(image, tempData, imgSize * sizeof(double));
		
	fftw_free(fftImageData);
	fftw_free(fftKernelData);
	fftw_free(tempData);
	
	lock.lock();
	fftw_destroy_plan(inToFPlan);
	fftw_destroy_plan(fToOutPlan);
	lock.unlock();
}

void FFTConvolver::Reverse(double* image, size_t imgWidth, size_t imgHeight)
{
	for(size_t y=0; y!=imgHeight/2; ++y)
	{
		size_t destY = imgHeight-1 - y;
		double* sourcePtr = &image[y*imgWidth];
		double* destPtr = &image[destY*imgWidth];
		for(size_t x=0; x!=imgWidth/2; ++x)
		{
			size_t destX = imgWidth-1 - x;
			std::swap(sourcePtr[x], destPtr[destX]);
		}
	}
}
