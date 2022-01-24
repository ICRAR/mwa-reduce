#ifndef FFT_CONVOLVER_H
#define FFT_CONVOLVER_H

#include <cstring>
#include <mutex>

class FFTConvolver {
 public:
  /**
   * Convolve an image with a smaller kernel. No preparation of either image is
   * needed.
   *
   * This function assumes that (n+1)/2 is the middle pixel for uneven image
   * sizes. In the case of even image sizes, the middle falls between two
   * pixels.
   */
  static void Convolve(class FFTWManager& fftw, float* image, size_t imgWidth,
                       size_t imgHeight, const float* kernel, size_t kernelSize,
                       size_t threadCount);

  /**
   * Convolve an image with a smaller kernel. No preparation of either image is
   * needed.
   *
   * This function assumes that n/2 is the middle pixel, and reverses the
   * kernel.
   */
  static void ReverseAndConvolve(class FFTWManager& fftw, float* image,
                                 size_t imgWidth, size_t imgHeight,
                                 const float* kernel, size_t kernelSize,
                                 size_t threadCount);

  /**
   * Prepare a smaller kernel for convolution with ConvolveSameSize. When the
   * kernel is used more often, it is more efficient to call PrepareKernel()
   * once and multiple times ConvolveSameSize(), than calling Convolve()
   * multiple times.
   */
  static void PrepareSmallKernel(float* dest, size_t imgWidth, size_t imgHeight,
                                 const float* kernel, size_t kernelSize,
                                 size_t threadCount);
  /**
   * Prepare a kernel for convolution with ConvolveSameSize(). The kernel should
   * be already of the same size as the image to be convolved, or otherwise
   * PrepareSmallKernel() should be used.
   */
  static void PrepareKernel(float* dest, const float* source, size_t imgWidth,
                            size_t imgHeight, size_t threadCount);

  /**
   * Convolve an image with an already prepared kernel of the same size.
   */
  static void ConvolveSameSize(class FFTWManager& fftw, float* image,
                               const float* kernel, size_t imgWidth,
                               size_t imgHeight, size_t threadCount);

  static void Reverse(float* image, size_t imgWidth, size_t imgHeight,
                      size_t threadCount);

 private:
};

#endif
