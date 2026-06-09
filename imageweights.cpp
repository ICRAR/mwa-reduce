#include "imageweights.h"

#include "msproviders/msprovider.h"

#include <aocommon/fits/fitswriter.h>
#include <aocommon/units/angle.h>
#include <aocommon/banddata.h>
#include <aocommon/multibanddata.h>

#include "wsclean/logger.h"

#include <cmath>
#include <iostream>
#include <cstring>

ImageWeights::ImageWeights(const WeightMode& weightMode, size_t imageWidth,
                           size_t imageHeight, double pixelScaleX,
                           double pixelScaleY, bool weightsAsTaper,
                           double superWeight)
    : _weightMode(weightMode),
      _imageWidth(round(double(imageWidth) / superWeight)),
      _imageHeight(round(double(imageHeight) / superWeight)),
      _pixelScaleX(pixelScaleX),
      _pixelScaleY(pixelScaleY),
      _totalSum(0.0),
      _isGriddingFinished(false),
      _weightsAsTaper(weightsAsTaper) {
  if (_imageWidth % 2 != 0) ++_imageWidth;
  if (_imageHeight % 2 != 0) ++_imageHeight;
  _grid.assign(_imageWidth * _imageHeight / 2, 0.0);
}

double ImageWeights::ApplyWeights(std::complex<float>* data, const bool* flags,
                                  double uTimesLambda, double vTimesLambda,
                                  size_t channelCount, double lowestFrequency,
                                  double frequencyStep) {
  double weightSum = 0.0;
  for (size_t ch = 0; ch != channelCount; ++ch) {
    if (flags[ch]) {
      data[ch] = 0.0;
    } else {
      double wavelength =
          frequencyToWavelength(lowestFrequency + frequencyStep * ch);
      double u = uTimesLambda / wavelength;
      double v = vTimesLambda / wavelength;
      double weight = GetWeight(u, v);
      weightSum += weight;
      data[ch] *= weight;
    }
  }
  return weightSum / channelCount;
}

void ImageWeights::Grid(casacore::MeasurementSet& ms,
                        const MSSelection& selection) {
  if (_isGriddingFinished)
    throw std::runtime_error("Grid() called after a call to FinishGridding()");
  const aocommon::MultiBandData bandData(ms.spectralWindow(),
                                         ms.dataDescription());
  casacore::ROScalarColumn<int> antenna1Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA1));
  casacore::ROScalarColumn<int> antenna2Column(
      ms, casacore::MS::columnName(casacore::MSMainEnums::ANTENNA2));
  casacore::ROScalarColumn<int> fieldIdColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::FIELD_ID));
  casacore::ROScalarColumn<double> timeColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
  casacore::ROArrayColumn<double> uvwColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::UVW));
  casacore::ROArrayColumn<float> weightColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::WEIGHT_SPECTRUM));
  casacore::ROArrayColumn<bool> flagColumn(
      ms, casacore::MS::columnName(casacore::MSMainEnums::FLAG));
  casacore::ROScalarColumn<int> dataDescIdColumn(
      ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));

  const casacore::IPosition shape(flagColumn.shape(0));

  bool hasWeights = false;
  if (weightColumn.isDefined(0) && !_weightsAsTaper) {
    casacore::IPosition weightShape = weightColumn.shape(0);
    hasWeights = (weightShape == shape);
  }

  const size_t polarizationCount = shape[0];

  casacore::Array<casacore::Complex> dataArr(shape);
  casacore::Array<bool> flagArr(shape);
  casacore::Array<float> weightArr(shape);
  size_t timestep = 0;
  double time = timeColumn(0);

  if (!hasWeights) weightArr.set(1.0);

  for (size_t row = 0; row != ms.nrow(); ++row) {
    const int a1 = antenna1Column(row), a2 = antenna2Column(row),
              fieldId = fieldIdColumn(row);
    if (time != timeColumn(row)) {
      ++timestep;
      time = timeColumn(row);
    }
    const casacore::Vector<double> uvw = uvwColumn(row);
    if (selection.IsSelected(fieldId, timestep, a1, a2, uvw)) {
      flagColumn.get(row, flagArr);
      if (hasWeights) weightColumn.get(row, weightArr);
      const aocommon::BandData& curBand = bandData[dataDescIdColumn(row)];

      bool* flagIter = flagArr.cbegin();
      float* weightIter = weightArr.cbegin();

      double uInM = uvw(0), vInM = uvw(1);
      if (vInM < 0.0) {
        uInM = -uInM;
        vInM = -vInM;
      }

      size_t startChannel, endChannel;
      if (selection.HasChannelRange()) {
        startChannel = selection.ChannelRangeStart();
        endChannel = selection.ChannelRangeEnd();
      } else {
        startChannel = 0;
        endChannel = curBand.ChannelCount();
      }

      for (size_t ch = startChannel; ch != endChannel; ++ch) {
        double u = uInM / curBand.ChannelWavelength(ch),
               v = vInM / curBand.ChannelWavelength(ch);
        int x, y;
        uvToXY(u, v, x, y);
        if (isWithinLimits(x, y)) {
          for (size_t p = 0; p != polarizationCount; ++p) {
            if (!*flagIter) {
              size_t index = (size_t)x + (size_t)y * _imageWidth;
              _grid[index] += *weightIter;
              _totalSum += *weightIter;
            }
            ++flagIter;
            ++weightIter;
          }
        } else {
          for (size_t p = 0; p != polarizationCount; ++p) {
            ++flagIter;
            ++weightIter;
          }
        }
      }
    }
  }
}

void ImageWeights::Grid(MSProvider& msProvider, const MSSelection& selection) {
  if (_isGriddingFinished)
    throw std::runtime_error("Grid() called after a call to FinishGridding()");
  size_t polarizationCount =
      (msProvider.Polarization() == aocommon::Polarization::Instrumental) ? 4
                                                                          : 1;
  if (_weightMode.RequiresGridding()) {
    SynchronizedMS ms(msProvider.MS());
    const aocommon::MultiBandData bandData(ms->spectralWindow(),
                                           ms->dataDescription());
    aocommon::MultiBandData selectedBand;
    if (selection.HasChannelRange())
      selectedBand = aocommon::MultiBandData(
          bandData, selection.ChannelRangeStart(), selection.ChannelRangeEnd());
    else
      selectedBand = bandData;
    aocommon::UVector<float> weightBuffer(selectedBand[0].ChannelCount() *
                                          polarizationCount);

    msProvider.Reset();
    while (msProvider.CurrentRowAvailable()) {
      double uInM, vInM, wInM;
      size_t dataDescId;
      msProvider.ReadMeta(uInM, vInM, wInM, dataDescId);
      msProvider.ReadWeights(weightBuffer.data());
      if (_weightsAsTaper) {
        for (float& w : weightBuffer) {
          if (w != 0.0) w = 1.0;
        }
      }
      const aocommon::BandData& curBand = selectedBand[dataDescId];
      if (vInM < 0.0) {
        uInM = -uInM;
        vInM = -vInM;
      }

      const float* weightIter = weightBuffer.data();
      for (size_t ch = 0; ch != curBand.ChannelCount(); ++ch) {
        double u = uInM / curBand.ChannelWavelength(ch),
               v = vInM / curBand.ChannelWavelength(ch);
        for (size_t p = 0; p != polarizationCount; ++p) {
          Grid(u, v, *weightIter);
          ++weightIter;
        }
      }

      msProvider.NextRow();
    }
  }
}

void ImageWeights::FinishGridding() {
  if (_isGriddingFinished)
    throw std::runtime_error("FinishGridding() called twice");
  _isGriddingFinished = true;

  switch (_weightMode.Mode()) {
    case WeightMode::BriggsWeighted: {
      double avgW = 0.0;
      for (aocommon::UVector<double>::const_iterator i = _grid.begin();
           i != _grid.end(); ++i)
        avgW += *i * *i;
      avgW /= _totalSum;
      double numeratorSqrt = 5.0 * exp10(-_weightMode.BriggsRobustness());
      double sSq = numeratorSqrt * numeratorSqrt / avgW;
      for (aocommon::UVector<double>::iterator i = _grid.begin();
           i != _grid.end(); ++i) {
        *i = 1.0 / (1.0 + *i * sSq);
      }
    } break;
    case WeightMode::UniformWeighted: {
      for (aocommon::UVector<double>::iterator i = _grid.begin();
           i != _grid.end(); ++i) {
        if (*i != 0.0)
          *i = 1.0 / *i;
        else
          *i = 0.0;
      }
    } break;
    case WeightMode::NaturalWeighted: {
      for (aocommon::UVector<double>::iterator i = _grid.begin();
           i != _grid.end(); ++i) {
        if (*i != 0.0) *i = 1.0;
      }
    } break;
  }
}

void ImageWeights::Grid(const std::complex<float>*, const bool* flags,
                        double uTimesLambda, double vTimesLambda,
                        size_t channelCount, double lowestFrequency,
                        double frequencyStep) {
  if (_isGriddingFinished)
    throw std::runtime_error("Grid() called after a call to FinishGridding()");
  for (size_t ch = 0; ch != channelCount; ++ch) {
    if (!flags[ch]) {
      if (vTimesLambda < 0.0) {
        uTimesLambda = -uTimesLambda;
        vTimesLambda = -vTimesLambda;
      }

      double wavelength =
          frequencyToWavelength(lowestFrequency + frequencyStep * ch);
      double x = round(uTimesLambda * _imageWidth * _pixelScaleX / wavelength +
                       _imageWidth / 2);
      double y = round(vTimesLambda * _imageHeight * _pixelScaleY / wavelength);
      if (x >= 0.0 && x < _imageWidth && y < _imageHeight / 2) {
        size_t index = (size_t)x + (size_t)y * _imageWidth;
        _grid[index] += 1.0;
      }
    }
  }
}

void ImageWeights::SetMinUVRange(double minUVInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  const double minSq = minUVInLambda * minUVInLambda;
  int halfWidth = _imageWidth / 2;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      int xi = int(x) - halfWidth;
      double u = double(xi) / (_imageWidth * _pixelScaleX);
      double v = double(y) / (_imageHeight * _pixelScaleY);
      if (u * u + v * v < minSq) *i = 0.0;
      ++i;
    }
  }
}

void ImageWeights::SetMaxUVRange(double maxUVInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  const double maxSq = maxUVInLambda * maxUVInLambda;
  int halfWidth = _imageWidth / 2;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      int xi = int(x) - halfWidth;
      double u = double(xi) / (_imageWidth * _pixelScaleX);
      double v = double(y) / (_imageHeight * _pixelScaleY);
      if (u * u + v * v > maxSq) *i = 0.0;
      ++i;
    }
  }
}

void ImageWeights::SetTukeyTaper(double transitionSizeInLambda,
                                 double maxUVInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  const double maxUVSq = maxUVInLambda * maxUVInLambda;
  const double transitionDistSq = (maxUVInLambda - transitionSizeInLambda) *
                                  (maxUVInLambda - transitionSizeInLambda);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double distSq = u * u + v * v;
      if (distSq > maxUVSq)
        *i = 0.0;
      else if (distSq > transitionDistSq) {
        *i *=
            tukeyFrom0ToN(maxUVInLambda - sqrt(distSq), transitionSizeInLambda);
      }
      ++i;
    }
  }
}

void ImageWeights::SetTukeyInnerTaper(double transitionSizeInLambda,
                                      double minUVInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  const double minUVSq = minUVInLambda * minUVInLambda;
  const double totalSizeSq = (minUVInLambda + transitionSizeInLambda) *
                             (minUVInLambda + transitionSizeInLambda);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double distSq = u * u + v * v;
      if (distSq < minUVSq)
        *i = 0.0;
      else if (distSq < totalSizeSq) {
        *i *=
            tukeyFrom0ToN(sqrt(distSq) - minUVInLambda, transitionSizeInLambda);
      }
      ++i;
    }
  }
}

void ImageWeights::SetEdgeTaper(double sizeInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  double maxU, maxV;
  xyToUV(_imageWidth, _imageHeight / 2, maxU, maxV);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      if (maxU - std::fabs(u) < sizeInLambda ||
          maxV - std::fabs(v) < sizeInLambda)
        *i = 0.0;
      ++i;
    }
  }
}

void ImageWeights::SetEdgeTukeyTaper(double transitionSizeInLambda,
                                     double edgeSizeInLambda) {
  aocommon::UVector<double>::iterator i = _grid.begin();
  double maxU, maxV;
  xyToUV(_imageWidth, _imageHeight / 2, maxU, maxV);
  double totalSize = transitionSizeInLambda + edgeSizeInLambda;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double u, v;
      xyToUV(x, y, u, v);
      double uDist = maxU - std::fabs(u);
      double vDist = maxV - std::fabs(v);
      if (uDist < edgeSizeInLambda || vDist < edgeSizeInLambda)
        *i = 0.0;
      else if (uDist < totalSize || vDist < totalSize) {
        double ru = uDist - edgeSizeInLambda;
        double rv = vDist - edgeSizeInLambda;
        if (ru > transitionSizeInLambda) ru = transitionSizeInLambda;
        if (rv > transitionSizeInLambda) rv = transitionSizeInLambda;
        *i *= tukeyFrom0ToN(ru, transitionSizeInLambda) *
              tukeyFrom0ToN(rv, transitionSizeInLambda);
      }
      ++i;
    }
  }
}

void ImageWeights::GetGrid(double* image) const {
  const double* srcPtr = _grid.data();
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    size_t yUpper = _imageHeight / 2 - 1 - y;
    size_t yLower = _imageHeight / 2 + y;
    double* upperRow = &image[yUpper * _imageWidth];
    double* lowerRow = &image[yLower * _imageWidth];
    for (size_t x = 0; x != _imageWidth; ++x) {
      upperRow[_imageWidth - x - 1] = *srcPtr;
      lowerRow[x] = *srcPtr;
      ++srcPtr;
    }
  }
}

void ImageWeights::Save(const string& filename) const {
  aocommon::UVector<double> image(_imageWidth * _imageHeight);
  GetGrid(&image[0]);
  aocommon::FitsWriter writer;
  writer.SetImageDimensions(_imageWidth, _imageHeight);
  writer.Write(filename, image.data());
}

void ImageWeights::RankFilter(double rankLimit, size_t windowSize) {
  aocommon::UVector<double> newGrid(_grid);
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double w = _grid[y * _imageWidth + x];
      if (w != 0.0) {
        double mean = windowMean(x, y, windowSize);
        if (w > mean * rankLimit)
          newGrid[y * _imageWidth + x] = mean * rankLimit;
      }
    }
  }
  _grid = newGrid;
}

void ImageWeights::SetGaussianTaper(double beamSize) {
  Logger::Debug << "Applying " << aocommon::units::Angle::ToNiceString(beamSize)
                << " Gaussian taper...\n";
  double halfPowerUV = 1.0 / (beamSize * 2.0 * M_PI);
  const long double sigmaToHP = 2.0L * sqrtl(2.0L * logl(2.0L));
  double minusTwoSigmaSq = halfPowerUV * sigmaToHP;
  Logger::Debug << "UV taper: " << minusTwoSigmaSq << '\n';
  minusTwoSigmaSq *= -2.0 * minusTwoSigmaSq;
  for (size_t y = 0; y != _imageHeight / 2; ++y) {
    for (size_t x = 0; x != _imageWidth; ++x) {
      double val = _grid[y * _imageWidth + x];
      if (val != 0.0) {
        double u, v;
        xyToUV(x, y, u, v);
        double gaus = exp((u * u + v * v) / minusTwoSigmaSq);
        _grid[y * _imageWidth + x] = val * gaus;
      }
    }
  }
}

double ImageWeights::windowMean(size_t x, size_t y, size_t windowSize) {
  size_t d = windowSize / 2;
  size_t x1, y1, x2, y2;
  if (x <= d)
    x1 = 0;
  else
    x1 = x - d;

  if (y <= d)
    y1 = 0;
  else
    y1 = y - d;

  if (x + d >= _imageWidth)
    x2 = _imageWidth;
  else
    x2 = x + d;

  if (y + d >= _imageHeight / 2)
    y2 = _imageHeight / 2;
  else
    y2 = y + d;

  size_t windowCount = 0;
  double windowSum = 0.0;
  for (size_t yi = y1; yi < y2; ++yi) {
    for (size_t xi = x1; xi < x2; ++xi) {
      double w = _grid[yi * _imageWidth + xi];
      if (w != 0.0) {
        ++windowCount;
        windowSum += w;
      }
    }
  }
  return windowSum / double(windowCount);
}
