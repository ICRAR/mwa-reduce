#ifndef IMAGE_WEIGHTS_H
#define IMAGE_WEIGHTS_H

#include <cstddef>
#include <complex>

#include <ms/MeasurementSets/MeasurementSet.h>

#include "uvector.h"

class ImageWeights
{
	public:
		ImageWeights(size_t imageWidth, size_t imageHeight, double pixelScale);
		
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
		double GetUniformWeight(double u, double v);
		double GetInverseTaperedWeight(double u, double v)
		{
			return sqrt(u*u + v*v);
		}

		void Grid(casa::MeasurementSet& ms);
		
		double ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

		void Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

	private:
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
		double _pixelScale;
		
		ao::uvector<double> _sum, _weight;
};

#endif
