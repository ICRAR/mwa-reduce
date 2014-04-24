#ifndef MULTI_SCALE_CLEAN_H
#define MULTI_SCALE_CLEAN_H

#include "cleanalgorithm.h"
#include "imageset.h"

#include "../uvector.h"

template<typename ImageSetType = clean_algorithms::PolarizedImageSet>
class MultiScaleClean : public CleanAlgorithm
{
public:
	void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain);
	
	static void MakeShapeFunction(double scaleSizeInPixels, ao::uvector<double>& output, size_t& n)
	{
		n = size_t(ceil(scaleSizeInPixels*0.5)*2.0)+1;
		output.resize(n * n);
		shapeFunction(n, output, scaleSizeInPixels);
	}
	
	static void Convolve(double* image, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize);
	static void PrepareKernel(double* dest, size_t imgWidth, size_t imgHeight, const double* kernel, size_t kernelSize);
	static void ConvolveSameSize(double* image, const double* kernel, size_t imgWidth, size_t imgHeight);
	
private:
	static void shapeFunction(size_t n, ao::uvector<double>& output2d, double scaleSizeInPixels)
	{
		if(scaleSizeInPixels == 0.0)
			output2d[0] = 1.0;
		else {
			double* outputPtr = output2d.data();
			for(int y=0; y!=int(n); ++y)
			{
				double dy = y - 0.5*(n-1);
				double dydy = dy * dy;
				for(int x=0; x!=int(n) ;++x)
				{
					double dx = x - 0.5*(n-1);
					double r = sqrt(dx*dx + dydy);
					*outputPtr = hannWindowFunction(r, n) * shapeFunction(r / scaleSizeInPixels);
					++outputPtr;
				}
			}
		}
	}
	
	static double hannWindowFunction(double x, size_t n)
	{
		return (x*2 <= n+1) ? (0.5 * (1.0 + cos(2.0*M_PI*x / double(n+1)))) : 0.0;
	}
	
	static double shapeFunction(double x)
	{
		return (x < 1.0) ? (1.0 - x*x) : 0.0;
	}
};

#endif
