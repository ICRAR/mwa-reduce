#ifndef IMAGE_WEIGHTS_H
#define IMAGE_WEIGHTS_H

#include <cstddef>
#include <complex>

#include <ms/MeasurementSets/MeasurementSet.h>

#include "uvector.h"
#include "inversionalgorithm.h"
#include "weightmode.h"
#include "msselection.h"

class ImageWeights
{
	public:
		ImageWeights(size_t imageWidth, size_t imageHeight, double pixelScaleX, double pixelScaleY, double superWeight=1.0);
		
		double GetWeight(double u, double v)
		{
			return GetUniformWeight(u ,v);
			//return GetCountWeight(u ,v);
			//return GetInverseTaperedWeight(u, v);
			//return GetNaturalWeight(u ,v);
		}
		double GetNaturalWeight(double u, double v) const
		{
			return 1.0;
		}
		double GetUniformWeight(double u, double v) const
		{
			double val = sumValue(u, v);
			if(val != 0.0)
				return 1.0 / val;
			else
				return 0.0;
		}
		double GetInverseTaperedWeight(double u, double v)
		{
			return sqrt(u*u + v*v);
		}
		double GetBriggsWeight(double u, double v) const
		{
			return sumValue(u, v);
		}

		void Grid(casa::MeasurementSet& ms, WeightMode weightMode, const MSSelection& selection);
		void Grid(class MSProvider& ms, WeightMode weightMode, const MSSelection& selection);
		
		double ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

		void Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

	private:
		ImageWeights(const ImageWeights&) { }
		void operator=(const ImageWeights&) { }
		
		double sumValue(double u, double v) const
		{
			if(v < 0.0) {
				u = -u;
				v = -v;
			}
			double x = round(u*_imageWidth*_pixelScaleX + _imageWidth/2);
			double y = round(v*_imageHeight*_pixelScaleY);
			if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
				return _sum[(size_t) x + (size_t) y*_imageWidth];
			else {
				return 0.0;
			}
		}
		
		template<typename T>
		static T frequencyToWavelength(const T frequency)
		{
			return speedOfLight() / frequency; 
		}
		static long double speedOfLight()
		{
			return 299792458.0L;
		}
		std::size_t _imageWidth, _imageHeight;
		double _pixelScaleX, _pixelScaleY;
		
		ao::uvector<double> _sum;
};

#endif
