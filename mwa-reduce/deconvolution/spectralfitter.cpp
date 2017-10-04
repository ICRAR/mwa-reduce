#include "spectralfitter.h"

#include "../polynomialfitter.h"
#include "../nlplfitter.h"

#include "../wsclean/logger.h"

#include <limits>

void SpectralFitter::FitAndEvaluate(double* values) const
{
	ao::uvector<double> terms;
	Fit(terms, values);
	Evaluate(values, terms);
}

void SpectralFitter::Fit(ao::uvector<double>& terms, const double* values) const
{
	switch(_mode)
	{
		default:
		case NoSpectralFitting:
			break;
			
		case PolynomialSpectralFitting: {
			PolynomialFitter fitter;
			double refFreq = ReferenceFrequency();
			for(size_t i=0; i!=_frequencies.size(); ++i)
			{
				fitter.AddDataPoint(_frequencies[i] / refFreq - 1.0, values[i], _weights[i]);
			}
			
			fitter.Fit(terms, _nTerms);
		} break;
		
		case LogPolynomialSpectralFitting: {
			NonLinearPowerLawFitter fitter;
			double refFreq = ReferenceFrequency();
			for(size_t i=0; i!=_frequencies.size(); ++i)
				fitter.AddDataPoint(_frequencies[i] / refFreq, values[i]);
			
			fitter.Fit(terms, _nTerms);
		} break;
	}
}

void SpectralFitter::Evaluate(double* values, const ao::uvector<double>& terms) const
{
	switch(_mode)
	{
		default:
		case NoSpectralFitting:
			break;
			
		case PolynomialSpectralFitting: {
			double refFreq = ReferenceFrequency();
			for(size_t i=0; i!=_frequencies.size(); ++i) {
				double newValue = PolynomialFitter::Evaluate(_frequencies[i] / refFreq - 1.0, terms);
				//std::cout << values[i] << "->" << newValue << ' ';
				values[i] = newValue;
			}
			//std::cout << '\n';
			
		} break;
		
		case LogPolynomialSpectralFitting: {
			double refFreq = ReferenceFrequency();
			for(size_t i=0; i!=_frequencies.size(); ++i) {
				double newValue = NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms, refFreq);
				//std::cout << values[i] << "->" << newValue << ' ';
				values[i] = newValue;
			}
			//std::cout << '\n';
			
		} break;
	}
}

double SpectralFitter::Evaluate(const ao::uvector<double>& terms, double frequency) const
{
	switch(_mode)
	{
		default:
		case NoSpectralFitting:
			throw std::runtime_error("Something is inconsistent: can't evaluate terms at frequency without fitting");
			
		case PolynomialSpectralFitting:
			return PolynomialFitter::Evaluate(frequency / ReferenceFrequency() - 1.0, terms);
		
		case LogPolynomialSpectralFitting:
			return NonLinearPowerLawFitter::Evaluate(frequency, terms, ReferenceFrequency());
	}
}
