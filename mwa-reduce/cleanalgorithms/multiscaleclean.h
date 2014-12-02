#ifndef MULTI_SCALE_CLEAN_H
#define MULTI_SCALE_CLEAN_H

#include "cleanalgorithm.h"
#include "imageset.h"
#include "simpleclean.h"

#include "../uvector.h"

template<typename ImageSetType>
class MultiScaleClean : public TypedCleanAlgorithm<ImageSetType>
{
public:
	MultiScaleClean(double beamSize, double pixelSizeX, double pixelSizeY) :
		_startScale(64.0 * beamSize / pixelSizeX),
		_minScale(0.0),
		_beamSize(beamSize),
		_pixelSizeX(pixelSizeX),
		_pixelSizeY(pixelSizeY)
	{ }
	
	virtual void ExecuteMajorIteration(ImageSetType& dataImage, ImageSetType& modelImage, std::vector<double*> psfImages, size_t width, size_t height, bool& reachedStopGain);
	
	static void MakeShapeFunction(double scaleSizeInPixels, ao::uvector<double>& output, size_t& n)
	{
		n = size_t(ceil(scaleSizeInPixels*0.5)*2.0)+1;
		output.resize(n * n);
		shapeFunction(n, output, scaleSizeInPixels);
	}
	
private:
	size_t _originalWidth, _originalHeight;
	size_t _rescaledWidth, _rescaledHeight;
	double
		_startScale, _minScale, _beamSize, _pixelSizeX, _pixelSizeY;
	ImageSetType *_dataImageLargeScale, *_dataImageNextScale, *_dataImageOriginal, *_modelImage;
	std::vector<double*> *_originalPsfs, *_scaledPsfs;
	
	void executeMajorIterationForScale(double currentScale, double nextScale, bool& reachedStopGain, bool& canCleanFurther);
	
	double scaleBiasFunction(double smallerScale, double largerScale) const
	{
		// From Cornwell 2008, "Multi-Scale CLEAN deconvolution of radio synthesis images"
		// S(alpha) = 1.0 - 0.6 alpha / alpha_Max
		return 1.0 - CleanAlgorithm::_multiscaleScaleBias * (smallerScale / largerScale);
	}
	
	static void shapeFunction(size_t n, ao::uvector<double>& output2d, double scaleSizeInPixels)
	{
		if(scaleSizeInPixels == 0.0)
			output2d[0] = 1.0;
		else {
			double sum = 0.0;
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
					sum += *outputPtr;
					++outputPtr;
				}
			}
			double normFactor = 1.0 / sum;
			for(ao::uvector<double>::iterator i=output2d.begin(); i!=output2d.end(); ++i)
				*i *= normFactor;
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
	
	struct CleanTask
	{
		size_t cleanCompX, cleanCompY;
		typename ImageSetType::Value peak;
	};
	
	struct CleanResult
	{
		CleanResult() : nextPeakX(0), nextPeakY(0), peakLevelUnnormalized(0.0)
		{ }
		size_t nextPeakX, nextPeakY;
		double peakLevelUnnormalized;
	};
	
	struct CleanThreadData
	{
		MultiScaleClean<ImageSetType>* parent;
		size_t startY, endY;
	};

	void findPeak(size_t& x, size_t& y, size_t startY, size_t stopY) const;
	
	void findPeak(size_t& x, size_t& y) const
	{
		findPeak(x, y, 0, _rescaledHeight);
	}
	
	std::string peakDescription(const ImageSetType& image, size_t x, size_t y, double rescaleFactor);
	
	void cleanThreadFunc(ao::lane<CleanTask>* taskLane, ao::lane<CleanResult>* resultLane, CleanThreadData cleanData);
	
	void subtractImage(double *image, const double *psf, size_t x, size_t y, double factor, size_t startY, size_t endY) const
	{
		SimpleClean::PartialSubtractImage(image, _rescaledWidth, _rescaledHeight, psf, _rescaledWidth, _rescaledHeight, x, y, factor, startY, endY);
	}
	
};

#endif
