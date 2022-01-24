#ifndef POLYNOMIAL_FITTER_H
#define POLYNOMIAL_FITTER_H

#include <aocommon/uvector.h>

#include <array>

class PolynomialFitter {
 public:
  typedef float num_t;
  void Clear() { _dataPoints.clear(); }

  void AddDataPoint(num_t x, num_t y, num_t w) {
    _dataPoints.emplace_back(std::array<num_t, 3>{{x, y, w}});
  }

  void Fit(aocommon::UVector<num_t>& terms, size_t nTerms);

  static num_t Evaluate(num_t x, const aocommon::UVector<num_t>& terms) {
    num_t val = terms[0];
    num_t f = 1.0;
    for (size_t i = 1; i != terms.size(); ++i) {
      f *= x;
      val += f * terms[i];
    }
    return val;
  }

  size_t size() const { return _dataPoints.size(); }

 private:
  aocommon::UVector<std::array<num_t, 3>> _dataPoints;
};

#endif
