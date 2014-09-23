#include "cleanalgorithm.h"

#include "../imagecoordinates.h"
#include "../modelsource.h"
#include "../model.h"

#include <unistd.h>

CleanAlgorithm::CleanAlgorithm() :
	_threshold(0.0),
	_subtractionGain(0.1),
	_stopGain(1.0),
	_cleanBorderRatio(0.05),
	_maxIter(500),
	_iterationNumber(0),
	_threadCount(sysconf(_SC_NPROCESSORS_ONLN)),
	_allowNegativeComponents(true),
	_stopOnNegativeComponent(false),
	_cleanAreas(0)
{
}

void CleanAlgorithm::GetModelFromImage(Model &model, const double* image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double phaseCentreDL, double phaseCentreDM, double spectralIndex, double refFreq, PolarizationEnum polarization)
{
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			double value = image[y*width + x];
			if(value != 0.0 && std::isfinite(value))
			{
				long double l, m;
				ImageCoordinates::XYToLM<long double>(x, y, pixelSizeX, pixelSizeY, width, height, l, m);
				l += phaseCentreDL; m += phaseCentreDM;
				ModelComponent component;
				long double ra, dec;
				ImageCoordinates::LMToRaDec<long double>(l, m, phaseCentreRA, phaseCentreDec, ra, dec);
				std::stringstream nameStr;
				nameStr << "component" << model.SourceCount();
				component.SetSED(SpectralEnergyDistribution(value, refFreq, spectralIndex, polarization));
				component.SetPosRA(ra);
				component.SetPosDec(dec);
				
				ModelSource source;
				source.SetName(nameStr.str());
				source.AddComponent(component);
				model.AddSource(source);
			}
		}
	}
}

void CleanAlgorithm::ResizeImage(double* dest, size_t newWidth, size_t newHeight, const double* source, size_t width, size_t height)
{
	size_t srcStartX = (width - newWidth) / 2, srcStartY = (height - newHeight) / 2;
	for(size_t y=0; y!=newHeight; ++y)
	{
		double* destPtr = dest + y * newWidth;
		const double* srcPtr = source + (y + srcStartY) * width + srcStartX;
		memcpy(destPtr, srcPtr, newWidth * sizeof(double));
	}
}

void CleanAlgorithm::RemoveNaNsInPSF(double* psf, size_t width, size_t height)
{
	double* endPtr = psf + width*height;
	while(psf != endPtr)
	{
		if(!std::isfinite(*psf)) *psf = 0.0;
		++psf;
	}
}

void CleanAlgorithm::CalculateFastCleanPSFSize(size_t& psfWidth, size_t& psfHeight, size_t imageWidth, size_t imageHeight)
{
	// With 2048 x 2048, the subtraction is already so quick that it is not really required to make the psf smaller
	if(imageWidth <= 2048)
		psfWidth = imageWidth;
	else if(imageWidth <= 4096)
		psfWidth = 2048;
	else
		psfWidth = imageWidth / 2;
	
	if(imageHeight <= 2048)
		psfHeight = imageHeight;
	else if(imageHeight <= 4096)
		psfHeight = 2048;
	else
		psfHeight = imageHeight / 2;
}
