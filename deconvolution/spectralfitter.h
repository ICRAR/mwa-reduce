#ifndef SPECTRAL_FITTER_H
#define SPECTRAL_FITTER_H

#include "../uvector.h"

enum SpectralFittingMode {
	NoSpectralFitting,
	PolynomialSpectralFitting,
	LogPolynomialSpectralFitting
};

class SpectralFitter
{
public:
	SpectralFitter(SpectralFittingMode mode, size_t nTerms) :
		_mode(mode), _nTerms(nTerms)
	{ }
	
	SpectralFittingMode Mode() const
	{
		return _mode;
	}
	
	void SetMode(SpectralFittingMode mode, size_t nTerms)
	{
		_mode = mode;
		_nTerms = nTerms;
	}
	
	void FitAndEvaluate(double* values) const;
	
	void Fit(ao::uvector<double>& terms, const double* values) const;
	
	void Evaluate(double* values, const ao::uvector<double>& terms) const;
	
	double Evaluate(const ao::uvector<double>& terms, double frequency) const;
	
	void SetFrequencies(const double* frequencies, const double* weights, size_t n)
	{
		_frequencies.assign(frequencies, frequencies+n);
		_weights.assign(weights, weights+n);
	}
	
	double Frequency(size_t index) const
	{
		return _frequencies[index];
	}
	
	double Weight(size_t index) const
	{
		return _weights[index];
	}
	
	size_t NTerms() const { return _nTerms; }
	
	size_t NFrequencies() const { return _frequencies.size(); }
	
	double ReferenceFrequency() const {
		return _frequencies[_frequencies.size()/2];
	}
	
private:
	enum SpectralFittingMode _mode;
	size_t _nTerms;
	ao::uvector<double> _frequencies, _weights;
};

#endif
