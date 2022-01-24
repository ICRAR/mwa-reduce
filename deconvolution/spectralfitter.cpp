#include "spectralfitter.h"

#include "../math/polynomialfitter.h"
#include "../math/nlplfitter.h"

#include "../wsclean/logger.h"

#include <limits>

namespace {
constexpr double kDefaultReferenceFrequency = 150e6;
}

void SpectralFitter::SetFrequencies(const double* frequencies,
                                    const num_t* weights, size_t n) {
  _frequencies.assign(frequencies, frequencies + n);
  _weights.assign(weights, weights + n);
  num_t weightSum = 0.0;
  _referenceFrequency = 0.0;
  for (size_t i = 0; i != n; ++i) {
    _referenceFrequency += _frequencies[i] * _weights[i];
    weightSum += _weights[i];
  }
  if (weightSum > 0.0)
    _referenceFrequency /= weightSum;
  else
    _referenceFrequency = kDefaultReferenceFrequency;
}

void SpectralFitter::Fit(aocommon::UVector<num_t>& terms, const num_t* values,
                         size_t x, size_t y) const {
  if (IsForced()) {
    forcedFit(terms, values, x, y);
  } else {
    switch (_mode) {
      default:
      case SpectralFittingMode::NoFitting:
        break;

      case SpectralFittingMode::Polynomial: {
        PolynomialFitter fitter;
        double refFreq = ReferenceFrequency();
        for (size_t i = 0; i != _frequencies.size(); ++i) {
          if (_weights[i] > 0.0) {
            fitter.AddDataPoint(_frequencies[i] / refFreq - 1.0, values[i],
                                _weights[i]);
          }
        }

        fitter.Fit(terms, _nTerms);
      } break;

      case SpectralFittingMode::LogPolynomial: {
        NonLinearPowerLawFitter fitter;
        double refFreq = ReferenceFrequency();
        for (size_t i = 0; i != _frequencies.size(); ++i) {
          if (_weights[i] > 0.0) {
            fitter.AddDataPoint(_frequencies[i] / refFreq, values[i]);
          }
        }

        fitter.Fit(terms, _nTerms);
      } break;
    }
  }
}

void SpectralFitter::forcedFit(aocommon::UVector<num_t>& terms,
                               const SpectralFitter::num_t* values, size_t x,
                               size_t y) const {
  terms.resize(_nTerms);
  // We need to find alpha such that
  // y[i] = A f(x[i], terms), with f the shape.
  // The least-squares fit is:
  // A = sum (y[i] w[i] f[i]) / sum (w[i] f[i]^2)
  // However, it turns out that finding the true least-squares solution for A
  // leads to unstable cleaning. This is because a LS constrained flux might
  // integrate to zero. If it does, the peak finding that uses integrated
  // values will again find the same peak (over and over...). Therefore,
  // we now use the linear average to estimate the flux:
  // A = sum (y[i] w[i]) / sum (w[i] f[i])
  // This is what is calculated below.
  terms[0] = 1.0;
  for (size_t term = 1; term != _nTerms; ++term) {
    const DImage& termImage = _forcedTerms[term - 1];
    terms[term] = termImage[x + y * termImage.Width()];
  }
  float aNumerator = 0.0;
  float aDivisor = 0.0;
  for (size_t i = 0; i != _frequencies.size(); ++i) {
    const float w = _weights[i];
    const float f = NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms,
                                                      ReferenceFrequency());
    if (w > 0.0) {
      aNumerator += w * values[i];
      aDivisor += w * f;
    }
  }
  terms[0] = aNumerator / aDivisor;
}

void SpectralFitter::Evaluate(num_t* values,
                              const aocommon::UVector<num_t>& terms) const {
  switch (_mode) {
    default:
    case SpectralFittingMode::NoFitting:
      break;

    case SpectralFittingMode::Polynomial: {
      const double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        values[i] =
            PolynomialFitter::Evaluate(_frequencies[i] / refFreq - 1.0, terms);
      }
    } break;

    case SpectralFittingMode::LogPolynomial: {
      const double refFreq = ReferenceFrequency();
      for (size_t i = 0; i != _frequencies.size(); ++i) {
        values[i] =
            NonLinearPowerLawFitter::Evaluate(_frequencies[i], terms, refFreq);
      }
    } break;
  }
}

SpectralFitter::num_t SpectralFitter::Evaluate(
    const aocommon::UVector<num_t>& terms, double frequency) const {
  switch (_mode) {
    default:
    case SpectralFittingMode::NoFitting:
      throw std::runtime_error(
          "Something is inconsistent: can't evaluate terms at frequency "
          "without fitting");

    case SpectralFittingMode::Polynomial:
      return PolynomialFitter::Evaluate(frequency / ReferenceFrequency() - 1.0,
                                        terms);

    case SpectralFittingMode::LogPolynomial:
      return NonLinearPowerLawFitter::Evaluate(frequency, terms,
                                               ReferenceFrequency());
  }
}
