#ifndef NLPL_FITTER_H
#define NLPL_FITTER_H

#include <cmath>
#include <memory>

#include "uvector.h"

/**
 * This class fits a power law to a set of points. Note that there is a
 * linear solution for this problem, but the linear solution requires
 * all values to be positive, which is not the case for e.g. spectral
 * energy distributions, because these have noise.
 * This fitter does not have this requirement.
 */
class NonLinearPowerLawFitter
{
public:
	NonLinearPowerLawFitter();
	
	~NonLinearPowerLawFitter();
	
	void AddDataPoint(double x, double y);
	
	void Fit(double& exponent, double& factor);
	
	void Fit(double& a, double& b, double& c);
	
	void Fit(ao::uvector<double>& terms, size_t nTerms);
	void FitStable(ao::uvector<double>& terms, size_t nTerms);
	
	void FastFit(double& exponent, double& factor);
	
	static double Evaluate(double x, const ao::uvector<double>& terms, double referenceFrequencyHz=1.0);
	
	static long double Evaluate(long double factor, long double exponent, long double frequencyHz)
	{
		return factor * powl(frequencyHz, exponent);
	}
		
	/*static double Term0ToFactor(double term0, double term1)
	{
		return exp10(term0); // + term1*log(NLPLFact));
	}
	
	static double FactorToTerm0(double factor, double term1)
	{
		return log10(factor); // - (term1*log(NLPLFact));
	}*/
	
private:
	void fit_implementation(ao::uvector<double>& terms, size_t nTerms);
	
	std::unique_ptr<class NLPLFitterData> _data;
};

#endif
