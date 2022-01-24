#ifndef IMAGE_WEIGHTS_H
#define IMAGE_WEIGHTS_H

#include <cstddef>
#include <complex>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <aocommon/uvector.h>
#include "inversionalgorithm.h"
#include "weightmode.h"
#include "msselection.h"

class ImageWeights {
 public:
  ImageWeights(const WeightMode& weightMode, size_t imageWidth,
               size_t imageHeight, double pixelScaleX, double pixelScaleY,
               bool weightsAsTaper, double superWeight);

  double GetWeight(double u, double v) const { return sampleGridValue(u, v); }

  void Grid(casacore::MeasurementSet& ms, const MSSelection& selection);
  void Grid(class MSProvider& ms, const MSSelection& selection);
  void Grid(double u, double v, double weight) {
    int x, y;
    uvToXY(u, v, x, y);

    if (isWithinLimits(x, y)) {
      size_t index = (size_t)x + (size_t)y * _imageWidth;
      _grid[index] += weight;
      _totalSum += weight;
    }
  }

  void FinishGridding();

  double ApplyWeights(std::complex<float>* data, const bool* flags,
                      double uTimesLambda, double vTimesLambda,
                      size_t channelCount, double lowestFrequency,
                      double frequencyStep);

  void Grid(const std::complex<float>* data, const bool* flags,
            double uTimesLambda, double vTimesLambda, size_t channelCount,
            double lowestFrequency, double frequencyStep);

  void SetMaxUVRange(double maxUVInLambda);
  void SetMinUVRange(double minUVInLambda);
  void SetGaussianTaper(double beamSize);
  void SetTukeyTaper(double transitionSizeInLambda, double maxUVInLambda);
  void SetTukeyInnerTaper(double transitionSizeInLambda, double minUVInLambda);
  void SetEdgeTaper(double sizeInLambda);
  void SetEdgeTukeyTaper(double transitionSizeInLambda,
                         double edgeSizeInLambda);

  void GetGrid(double* image) const;
  void Save(const std::string& filename) const;

  void RankFilter(double rankLimit, size_t windowSize);

  size_t Width() const { return _imageWidth; }
  size_t Height() const { return _imageHeight; }

 private:
  ImageWeights(const ImageWeights&) = delete;
  void operator=(const ImageWeights&) = delete;

  void uvToXY(double u, double v, int& x, int& y) const {
    if (v < 0.0) {
      u = -u;
      v = -v;
    }
    x = int(floor(u * _imageWidth * _pixelScaleX + _imageWidth / 2));
    y = int(floor(v * _imageHeight * _pixelScaleY));
  }

  void xyToUV(int x, int y, double& u, double& v) const {
    if (y < 0.0) {
      x = -x;
      y = -y;
    }

    u = double(x - int(_imageWidth / 2)) / double(_imageWidth * _pixelScaleX);
    v = double(y) / double(_imageHeight * _pixelScaleY);
  }

  bool isWithinLimits(int x, int y) const {
    return x >= 0 && x < int(_imageWidth) && y < int(_imageHeight / 2);
  }

  double sampleGridValue(double u, double v) const {
    int x, y;
    uvToXY(u, v, x, y);
    if (isWithinLimits(x, y))
      return _grid[(size_t)x + (size_t)y * _imageWidth];
    else {
      return 0.0;
    }
  }

  double windowMean(size_t x, size_t y, size_t windowSize);

  /**
   * Returns Tukey tapering function. This function is
   * 0 when x=0 and 1 when x=n.
   */
  double tukeyFrom0ToN(double x, double n) {
    return 0.5 * (1.0 + cos((M_PI / n) * (x - n)));
  }

  template <typename T>
  static T frequencyToWavelength(const T frequency) {
    return speedOfLight() / frequency;
  }
  static long double speedOfLight() { return 299792458.0L; }
  const WeightMode _weightMode;
  std::size_t _imageWidth, _imageHeight;
  const double _pixelScaleX, _pixelScaleY;

  aocommon::UVector<double> _grid;
  double _totalSum;
  bool _isGriddingFinished, _weightsAsTaper;
};

#endif
