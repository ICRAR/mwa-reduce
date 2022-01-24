#include "fftconvolver.h"
#include "fftkernels.h"

#include "fftwmanager.h"
#include <aocommon/uvector.h>
#include <aocommon/staticfor.h>

#include <fftw3.h>

#include <complex>
#include <stdexcept>

void FFTConvolver::Convolve(FFTWManager& fftw, float* image, size_t imgWidth,
                            size_t imgHeight, const float* kernel,
                            size_t kernelSize, size_t threadCount) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);
  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize, threadCount);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight,
                   threadCount);
}

void FFTConvolver::ReverseAndConvolve(class FFTWManager& fftw, float* image,
                                      size_t imgWidth, size_t imgHeight,
                                      const float* kernel, size_t kernelSize,
                                      size_t threadCount) {
  aocommon::UVector<float> scaledKernel(imgWidth * imgHeight, 0.0);

  PrepareSmallKernel(scaledKernel.data(), imgWidth, imgHeight, kernel,
                     kernelSize, threadCount);
  ConvolveSameSize(fftw, image, scaledKernel.data(), imgWidth, imgHeight,
                   threadCount);
}

void FFTConvolver::PrepareSmallKernel(float* dest, size_t imgWidth,
                                      size_t imgHeight, const float* kernel,
                                      size_t kernelSize, size_t threadCount) {
  if (kernelSize > imgWidth || kernelSize > imgHeight)
    throw std::runtime_error("Kernel size > image dimension");
  aocommon::StaticFor<size_t> loop(threadCount);
  loop.Run(0, kernelSize / 2, [&](size_t yStart, size_t yEnd) {
    const float* kernelIter = &kernel[yStart * kernelSize];
    for (size_t y = yStart; y != yEnd; ++y) {
      size_t destY = imgHeight - kernelSize / 2 + y;
      size_t firstX = imgWidth - kernelSize / 2;
      float* destIter = &dest[destY * imgWidth + firstX];
      for (size_t x = 0; x != kernelSize / 2; ++x) {
        *destIter = *kernelIter;
        ++kernelIter;
        ++destIter;
      }
      destIter = &dest[destY * imgWidth];
      for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
        *destIter = *kernelIter;
        ++kernelIter;
        ++destIter;
      }
    }
  });
  loop.Run(kernelSize / 2, kernelSize, [&](size_t yStart, size_t yEnd) {
    const float* kernelIter = &kernel[yStart * kernelSize];
    for (size_t y = yStart; y != yEnd; ++y) {
      size_t firstX = imgWidth - kernelSize / 2;
      float* destIter = &dest[firstX + (y - kernelSize / 2) * imgWidth];
      for (size_t x = 0; x != kernelSize / 2; ++x) {
        *destIter = *kernelIter;
        ++kernelIter;
        ++destIter;
      }
      destIter = &dest[(y - kernelSize / 2) * imgWidth];
      for (size_t x = kernelSize / 2; x != kernelSize; ++x) {
        *destIter = *kernelIter;
        ++kernelIter;
        ++destIter;
      }
    }
  });
}

void FFTConvolver::PrepareKernel(float* dest, const float* source,
                                 size_t imgWidth, size_t imgHeight,
                                 size_t threadCount) {
  aocommon::StaticFor<size_t> loop(threadCount);
  loop.Run(0, imgHeight / 2, [&](size_t yStart, size_t yEnd) {
    const float* sourceIter = &source[yStart * imgWidth];
    for (size_t y = yStart; y != yEnd; ++y) {
      size_t destY = imgHeight - imgHeight / 2 + y;
      size_t firstX = imgWidth - imgWidth / 2;
      float* destIter = &dest[destY * imgWidth + firstX];
      for (size_t x = 0; x != imgWidth / 2; ++x) {
        *destIter = *sourceIter;
        ++sourceIter;
        ++destIter;
      }
      destIter = &dest[destY * imgWidth];
      for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
        *destIter = *sourceIter;
        ++sourceIter;
        ++destIter;
      }
    }
  });
  loop.Run(imgHeight / 2, imgHeight, [&](size_t yStart, size_t yEnd) {
    const float* sourceIter = &source[yStart * imgWidth];
    for (size_t y = yStart; y != yEnd; ++y) {
      size_t firstX = imgWidth - imgWidth / 2;
      float* destIter = &dest[firstX + (y - imgHeight / 2) * imgWidth];
      for (size_t x = 0; x != imgWidth / 2; ++x) {
        *destIter = *sourceIter;
        ++sourceIter;
        ++destIter;
      }
      destIter = &dest[(y - imgHeight / 2) * imgWidth];
      for (size_t x = imgWidth / 2; x != imgWidth; ++x) {
        *destIter = *sourceIter;
        ++sourceIter;
        ++destIter;
      }
    }
  });
}

void FFTConvolver::ConvolveSameSize(FFTWManager& fftw, float* image,
                                    const float* kernel, size_t imgWidth,
                                    size_t imgHeight, size_t threadCount) {
  const size_t imgSize = imgWidth * imgHeight;
  const size_t complexWidth = imgWidth / 2 + 1;
  const size_t complexSize = complexWidth * imgHeight;
  float* tempData = fftwf_alloc_real(imgSize);
  fftwf_complex* fftImageData = fftwf_alloc_complex(complexSize);
  fftwf_complex* fftKernelData = fftwf_alloc_complex(complexSize);

  std::unique_lock<std::mutex> lock(fftw.Mutex());
  fftwf_plan plan_r2c =
      fftwf_plan_dft_r2c_1d(imgWidth, nullptr, nullptr, FFTW_ESTIMATE);
  fftwf_plan plan_c2c_forward = fftwf_plan_dft_1d(imgHeight, nullptr, nullptr,
                                                  FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_plan plan_c2c_backward = fftwf_plan_dft_1d(
      imgHeight, nullptr, nullptr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_plan plan_c2r =
      fftwf_plan_dft_c2r_1d(imgWidth, nullptr, nullptr, FFTW_ESTIMATE);
  lock.unlock();

  aocommon::StaticFor<size_t> loop(threadCount);

  fft2f_r2c_composite(plan_r2c, plan_c2c_forward, imgHeight, imgWidth, image,
                      fftImageData, loop);

  std::copy_n(kernel, imgSize, tempData);
  fft2f_r2c_composite(plan_r2c, plan_c2c_forward, imgHeight, imgWidth, tempData,
                      fftKernelData, loop);

  float fact = 1.0 / imgSize;
  loop.Run(0, imgHeight, [&](size_t yStart, size_t yEnd) {
    for (size_t y = yStart; y != yEnd; ++y) {
      for (size_t x = 0; x != complexWidth; ++x) {
        size_t i = y * complexWidth + x;
        reinterpret_cast<std::complex<float>*>(fftImageData)[i] *=
            fact * reinterpret_cast<std::complex<float>*>(fftKernelData)[i];
      }
    }
  });

  fft2f_c2r_composite(plan_c2c_backward, plan_c2r, imgHeight, imgWidth,
                      fftImageData, image, loop);

  fftwf_free(fftImageData);
  fftwf_free(fftKernelData);
  fftwf_free(tempData);

  lock.lock();
  fftwf_destroy_plan(plan_r2c);
  fftwf_destroy_plan(plan_c2c_forward);
  fftwf_destroy_plan(plan_c2c_backward);
  fftwf_destroy_plan(plan_c2r);
  lock.unlock();
}

void FFTConvolver::Reverse(float* image, size_t imgWidth, size_t imgHeight,
                           size_t threadCount) {
  aocommon::StaticFor<size_t> loop(threadCount);
  loop.Run(0, imgHeight / 2, [&](size_t yStart, size_t yEnd) {
    for (size_t y = yStart; y != yEnd; ++y) {
      size_t destY = imgHeight - 1 - y;
      float* sourcePtr = &image[y * imgWidth];
      float* destPtr = &image[destY * imgWidth];
      for (size_t x = 0; x != imgWidth / 2; ++x) {
        size_t destX = imgWidth - 1 - x;
        std::swap(sourcePtr[x], destPtr[destX]);
      }
    }
  });
}
