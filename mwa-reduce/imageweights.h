#ifndef IMAGE_WEIGHTS_H
#define IMAGE_WEIGHTS_H

#include <cstddef>
#include <complex>

#include <ms/MeasurementSets/MeasurementSet.h>

#include "uvector.h"

class ImageWeights
{
	public:
		ImageWeights(size_t imageWidth, size_t imageHeight, size_t channelCount, double pixelScale, double lowestFrequency, double frequencyStep);
		
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
		
		double ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda);

		void Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda);

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
		std::size_t _imageWidth, _imageHeight, _channelCount;
		double _pixelScale, _lowestFrequency, _frequencyStep;
		
		uvector<double> _sum, _weight;
};

#endif
