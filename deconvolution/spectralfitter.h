#ifndef SPECTRAL_FITTER_H
#define SPECTRAL_FITTER_H

#include <aocommon/uvector.h>
#include <vector>

#include "../image.h"

enum class SpectralFittingMode { NoFitting, Polynomial, LogPolynomial };

class SpectralFitter {
 public:
  typedef float num_t;

  SpectralFitter(SpectralFittingMode mode, size_t nTerms)
      : _mode(mode), _nTerms(nTerms) {}

  SpectralFittingMode Mode() const { return _mode; }

  void SetMode(SpectralFittingMode mode, size_t nTerms) {
    _mode = mode;
    _nTerms = nTerms;
  }

  /**
   * Fit an array of values to a curve.
   *
   * The type of the curve is set in the constructor or with @ref SetMode().
   * The coordinates are used in case the forced term fitting mode is used, in
   * which case it is used to look up the spectral index (or other terms) from
   * a specified image.
   *
   * @param [out] terms will hold the fitted terms. The meaning of these terms
   * depends on the fitted curve type, and are relative to the reference
   * frequency.
   * @param values array of size @ref NFrequencies() with the values to be
   * fitted. values[i] should correspond Frequency(i) and Weight(i).
   * @param x a pixel index giving the horizontal position
   * @param y a pixel index giving the vertical position
   */
  void Fit(aocommon::UVector<num_t>& terms, const num_t* values, size_t x,
           size_t y) const;

  /**
   * Evaluate the curve at the initialized frequencies.
   *
   * @param values array of size @ref NFrequencies() that will be filled with
   * curve values.
   * @param terms array of size @ref NTerms() with previously fitted terms.
   */
  void Evaluate(num_t* values, const aocommon::UVector<num_t>& terms) const;

  /**
   * Evaluate the curve at a specified frequency.
   *
   * @param terms array of size @ref NTerms() with previously fitted terms.
   * @param frequency Frequency in Hz.
   */
  num_t Evaluate(const aocommon::UVector<num_t>& terms, double frequency) const;

  /**
   * Fit an array of values to a curve, and replace those values
   * with the curve values. This function combines @ref Fit()
   * and @ref Evaluate().
   *
   * @param terms is a UVector of any size, that is used to store the terms.
   * Having this parameter explicitly is useful to avoid repeated allocation,
   * to temporarily store the terms: This function is used in the reasonably
   * critical loops inside deconvolution. It will be resized to @ref NTerms().
   */
  void FitAndEvaluate(num_t* values, size_t x, size_t y,
                      aocommon::UVector<num_t>& terms) const {
    Fit(terms, values, x, y);
    Evaluate(values, terms);
  }

  void SetFrequencies(const double* frequencies, const num_t* weights,
                      size_t n);

  double Frequency(size_t index) const { return _frequencies[index]; }

  num_t Weight(size_t index) const { return _weights[index]; }

  size_t NTerms() const { return _nTerms; }

  size_t NFrequencies() const { return _frequencies.size(); }

  double ReferenceFrequency() const { return _referenceFrequency; }

  void SetForcedImages(std::vector<DImage>&& images) {
    _forcedTerms = std::move(images);
  }

  bool IsForced() const { return !_forcedTerms.empty(); }

 private:
  void forcedFit(aocommon::UVector<num_t>& terms, const num_t* values, size_t x,
                 size_t y) const;

  enum SpectralFittingMode _mode;
  size_t _nTerms;
  aocommon::UVector<double> _frequencies;
  aocommon::UVector<num_t> _weights;
  double _referenceFrequency;
  std::vector<DImage> _forcedTerms;
};

#endif
