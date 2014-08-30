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
		ImageWeights(const WeightMode& weightMode, size_t imageWidth, size_t imageHeight, double pixelScaleX, double pixelScaleY, double superWeight=1.0);
		
		double GetNaturalWeight(double u, double v) const
		{
			return 1.0;
		}
		
		double GetUniformWeight(double u, double v) const
		{
			double val = sampleGridValue(u, v);
			if(val != 0.0)
				return 1.0 / val;
			else
				return 0.0;
		}
		
		double GetInverseTaperedWeight(double u, double v) const
		{
			return sqrt(u*u + v*v);
		}
		
		double GetWeight(double u, double v) const
		{
			switch(_weightMode.Mode())
			{
				case WeightMode::UniformWeighted:
				case WeightMode::BriggsWeighted: return sampleGridValue(u, v);
				case WeightMode::DistanceWeighted: return GetInverseTaperedWeight(u, v);
				default:
				case WeightMode::NaturalWeighted: return GetNaturalWeight(u, v);
			}
		}

		void Grid(casa::MeasurementSet& ms, const MSSelection& selection);
		void Grid(class MSProvider& ms, const MSSelection& selection);
		
		void FinishGridding();
		
		double ApplyWeights(std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

		void Grid(const std::complex<float> *data, const bool *flags, double uTimesLambda, double vTimesLambda, size_t channelCount, double lowestFrequency, double frequencyStep);

	private:
		ImageWeights(const ImageWeights&) :
			_weightMode(WeightMode::NaturalWeighted),
			_imageWidth(0),
			_imageHeight(0),
			_pixelScaleX(0.0),
			_pixelScaleY(0.0),
			_totalSum(0.0)
		{ }
		void operator=(const ImageWeights&) { }
		
		double sampleGridValue(double u, double v) const
		{
			if(v < 0.0) {
				u = -u;
				v = -v;
			}
			double x = round(u*_imageWidth*_pixelScaleX + _imageWidth/2);
			double y = round(v*_imageHeight*_pixelScaleY);
			if(x >= 0.0 && x < _imageWidth && y < _imageHeight/2)
				return _grid[(size_t) x + (size_t) y*_imageWidth];
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
		const WeightMode _weightMode;
		std::size_t _imageWidth, _imageHeight;
		const double _pixelScaleX, _pixelScaleY;
		
		ao::uvector<double> _grid;
		double _totalSum;
		bool _isGriddingFinished;
};

#endif
