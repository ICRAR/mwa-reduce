#ifndef NLPL_FITTER_H
#define NLPL_FITTER_H

#include <cmath>
#include <memory>

#include <aocommon/uvector.h>

/**
 * This class fits a power law to a set of points. Note that there is a
 * linear solution for this problem, but the linear solution requires
 * all values to be positive, which is not the case for e.g. spectral
 * energy distributions, because these have noise.
 * This fitter does not have this requirement.
 */
class NonLinearPowerLawFitter {
 public:
  typedef float num_t;

  NonLinearPowerLawFitter();

  ~NonLinearPowerLawFitter();

  void AddDataPoint(num_t x, num_t y);

  void Fit(num_t& exponent, num_t& factor);

  void Fit(num_t& a, num_t& b, num_t& c);

  void Fit(aocommon::UVector<num_t>& terms, size_t nTerms);
  void FitStable(aocommon::UVector<num_t>& terms, size_t nTerms);

  void FastFit(num_t& exponent, num_t& factor);

  static num_t Evaluate(num_t x, const aocommon::UVector<num_t>& terms,
                        num_t referenceFrequencyHz = 1.0);

  static long double Evaluate(long double factor, long double exponent,
                              long double frequencyHz) {
    return factor * powl(frequencyHz, exponent);
  }

 private:
  void fit_implementation(aocommon::UVector<num_t>& terms, size_t nTerms);

  std::unique_ptr<class NLPLFitterData> _data;
};

#endif
