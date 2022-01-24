#ifndef MODEL_MEASUREMENT_H
#define MODEL_MEASUREMENT_H

//#define EXTRA_ASSERTIONS 1

#include <cstring>

#include <aocommon/polarization.h>

class Measurement {
 public:
  Measurement() : _frequencyHz(0.0), _bandWidthHz(0.0) {
    for (size_t p = 0; p != 4; ++p) {
      _fluxDensities[p] = 0;
      _fluxDensityStddevs[p] = 0;
    }
  }

  Measurement(const Measurement &source)
      : _frequencyHz(source._frequencyHz), _bandWidthHz(source._bandWidthHz) {
    memcpy(_fluxDensities, source._fluxDensities, sizeof(long double) * 4);
    memcpy(_fluxDensityStddevs, source._fluxDensityStddevs,
           sizeof(long double) * 4);
  }

  Measurement &operator=(const Measurement &source) {
    _frequencyHz = source._frequencyHz;
    _bandWidthHz = source._bandWidthHz;
    memcpy(_fluxDensities, source._fluxDensities, sizeof(long double) * 4);
    memcpy(_fluxDensityStddevs, source._fluxDensityStddevs,
           sizeof(long double) * 4);
    return *this;
  }

  void operator+=(const Measurement &rhs) {
    for (size_t p = 0; p != 4; ++p) {
      _fluxDensities[p] += rhs._fluxDensities[p];
    }
  }

  void AverageWidth(const Measurement &rhs) {
    for (size_t p = 0; p != 4; ++p) {
      _fluxDensities[p] = (_fluxDensities[p] + rhs._fluxDensities[p]) * 0.5;
    }
  }

  void AverageWidth(const Measurement &rhs, double weight) {
    for (size_t p = 0; p != 4; ++p) {
      _fluxDensities[p] =
          (_fluxDensities[p] * (1.0 - weight) + rhs._fluxDensities[p] * weight);
    }
  }

  long double FrequencyHz() const { return _frequencyHz; }

  void SetFrequencyHz(long double frequencyHz) { _frequencyHz = frequencyHz; }

  long double FluxDensity(aocommon::PolarizationEnum polarization) const {
#ifdef EXTRA_ASSERTIONS
    if (!aocommon::Polarization::IsStokes(polarization))
      throw std::runtime_error("Cannot store specified polarization in model");
#endif
    return _fluxDensities[aocommon::Polarization::StokesToIndex(polarization)];
  }

  long double FluxDensityFromIndex(size_t polarizationIndex) const {
    return _fluxDensities[polarizationIndex];
  }

  void SetFluxDensityFromIndex(size_t polarizationIndex, long double flux) {
    _fluxDensities[polarizationIndex] = flux;
  }

  void SetFluxDensity(aocommon::PolarizationEnum polarization,
                      long double flux) {
    _fluxDensities[aocommon::Polarization::StokesToIndex(polarization)] = flux;
  }

  void SetZeroExceptSinglePol(aocommon::PolarizationEnum polarization,
                              long double flux) {
    _fluxDensities[0] = 0.0;
    _fluxDensities[1] = 0.0;
    _fluxDensities[2] = 0.0;
    _fluxDensities[3] = 0.0;
    _fluxDensities[aocommon::Polarization::StokesToIndex(polarization)] = flux;
#ifdef EXTRA_ASSERTIONS
    if (!aocommon::Polarization::IsStokes(polarization))
      throw std::runtime_error("Cannot store specified polarization in model");
#endif
  }

  void SetFluxDensityStddevFromIndex(size_t polarizationIndex,
                                     long double stddev) {
    _fluxDensityStddevs[polarizationIndex] = stddev;
  }

  void SetFluxDensityStddev(aocommon::PolarizationEnum polarization,
                            long double stddev) {
    _fluxDensityStddevs[aocommon::Polarization::StokesToIndex(polarization)] =
        stddev;
#ifdef EXTRA_ASSERTIONS
    if (!aocommon::Polarization::IsStokes(polarization))
      throw std::runtime_error("Cannot store specified polarization in model");
#endif
  }

  void SetBandWidthHz(double bandwidthHz) { _bandWidthHz = bandwidthHz; }

  void ToStream(std::ostream &s) const {
    s << "    measurement {\n"
         "      frequency "
      << (_frequencyHz / 1000000.0)
      << " MHz\n"
         "      fluxdensity Jy "
      << _fluxDensities[0] << ' ' << _fluxDensities[1] << ' '
      << _fluxDensities[2] << ' ' << _fluxDensities[3] << '\n';
    if (_bandWidthHz > 0.0) s << "      bandwidth " << _bandWidthHz << " Hz\n";
    s << "    }\n";
  }

 private:
  double _frequencyHz;
  double _bandWidthHz;
  long double _fluxDensities[4];
  long double _fluxDensityStddevs[4];
};

#endif
