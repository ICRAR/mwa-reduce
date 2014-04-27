#ifndef FFT_CONVOLVER_H
#define FFT_CONVOLVER_H

#include <cstring>

class FFTConvolver {
public:
	static void Convolve(double* image, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize);
	static void PrepareKernel(double* dest, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize);
	static void PrepareForConvolution(double* dest, const double* source, size_t imgWidth, size_t imgHeight);
	static void ConvolveSameSize(double* image, const double* kernel, size_t imgWidth, size_t imgHeight);
};

#endif
