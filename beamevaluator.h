#ifndef BEAM_EVALUATOR_H
#define BEAM_EVALUATOR_H

// 202008 : GPU version of beam2016 model added (developed by Maciej Cytowski
// and Di Pietrantonio from PAWSEY)

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>

#include <memory>

#include "beam/tilebeam.h"

#include <aocommon/matrix2x2.h>

#include "model/modelsource.h"
#include "model/model.h"

class BeamEvaluator {
 public:
  typedef TileBeam::PrecalcPosInfo PrecalcPosInfo;

  BeamEvaluator(casacore::MeasurementSet& ms, bool reportDelays,
                const std::string& mwaPath);

  void EvaluateApparentToAbsGain(double ra, double dec,
                                 std::complex<double>* gains) {
    EvaluateApparentToAbsGain(ra, dec, _frequency, gains);
  }

  void EvaluateAbsToApparentGain(double ra, double dec,
                                 std::complex<double>* gains) {
    EvaluateAbsToApparentGain(ra, dec, _frequency, gains);
  }

  void EvaluateApparentToAbsGain(double ra, double dec, double frequency,
                                 std::complex<double>* gains) {
    EvaluateAbsToApparentGain(ra, dec, frequency, gains);
    aocommon::Matrix2x2::Invert(gains);
  }

  void EvaluateAbsToApparentGain(double ra, double dec, double frequency,
                                 std::complex<double>* gains);

  template <typename NumType>
  void AbsToApparent(double ra, double dec,
                     std::complex<NumType>* pixelValues) {
    AbsToApparent<NumType>(ra, dec, _frequency, pixelValues);
  }

  template <typename NumType>
  void AbsToApparent(double ra, double dec, double frequency,
                     std::complex<NumType>* pixelValues) {
    std::complex<NumType> gains[4], temp[4];
    EvaluateAbsToApparentGain(ra, dec, frequency, gains);

    // Calculate A D A^H
    aocommon::Matrix2x2::ATimesB(temp, gains, pixelValues);
    aocommon::Matrix2x2::ATimesHermB(pixelValues, temp, gains);
  }

  template <typename NumType>
  void AbsToApparent(double ra, double dec, double frequency,
                     NumType* pixelValues) {
    std::complex<NumType> input[4], gains[4], temp[4];
    EvaluateAbsToApparentGain(ra, dec, frequency, gains);
    input[0] = pixelValues[0];
    input[1] = pixelValues[1];
    input[2] = pixelValues[2];
    input[3] = pixelValues[3];
    // Calculate A D A^H
    aocommon::Matrix2x2::ATimesB(temp, gains, input);
    aocommon::Matrix2x2::ATimesHermB(input, temp, gains);
    pixelValues[0] = input[0].real();
    pixelValues[1] = input[1].real();
    pixelValues[2] = input[2].real();
    pixelValues[3] = input[3].real();
  }

  template <typename NumType>
  void ApparentToAbs(double ra, double dec, std::complex<NumType>* data) {
    ApparentToAbs<NumType>(ra, dec, _frequency, data);
  }

  template <typename NumType>
  void ApparentToAbs(double ra, double dec, double frequency,
                     std::complex<NumType>* data) {
    std::complex<double> gains[4], temp[4];
    EvaluateApparentToAbsGain(ra, dec, frequency, gains);

    aocommon::Matrix2x2::ATimesB(temp, gains, data);
    aocommon::Matrix2x2::ATimesHermB(data, temp, gains);
  }

  template <typename NumType>
  void ApparentToAbs(double ra, double dec, NumType* pixelValues) {
    ApparentToAbs<NumType>(ra, dec, _frequency, pixelValues);
  }

  template <typename NumType>
  void ApparentToAbs(double ra, double dec, double frequency,
                     NumType* pixelValues) {
    std::complex<double> input[4], gains[4], temp[4];
    EvaluateApparentToAbsGain(ra, dec, frequency, gains);
    input[0] = pixelValues[0];
    input[1] = pixelValues[1];
    input[2] = pixelValues[2];
    input[3] = pixelValues[3];
    // Calculate A D A^H
    aocommon::Matrix2x2::ATimesB(temp, gains, input);
    aocommon::Matrix2x2::ATimesHermB(input, temp, gains);
    pixelValues[0] = input[0].real();
    pixelValues[1] = input[1].real();
    pixelValues[2] = input[2].real();
    pixelValues[3] = input[3].real();
  }

  void SetTime(const casacore::MEpoch& time) { _time = time; }
  const casacore::MEpoch& Time() { return _time; }

  void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, double raRad,
                                double decRad) {
    _tileBeam->PrecalculatePositionInfo(posInfo, _time, _ant1Pos, raRad,
                                        decRad);
  }

#ifndef CUDA_SUPPORT
  void EvaluateApparentToAbsGain(const PrecalcPosInfo& posInfo,
                                 double frequency,
                                 std::complex<double>* gains) {
    EvaluateAbsToApparentGain(posInfo, frequency, gains);
    aocommon::Matrix2x2::Invert(gains);
  }

  void EvaluateAbsToApparentGain(const PrecalcPosInfo& posInfo,
                                 double frequency,
                                 std::complex<double>* gains) {
    _tileBeam->ArrayResponse(posInfo, frequency, gains);
  }
#else
  void EvaluateAbsToApparentGain(const std::vector<PrecalcPosInfo>& posInfo,
                                 std::complex<double>* gains,
                                 size_t channelStart, size_t endChannel,
                                 size_t channelCount, size_t startFrequency,
                                 size_t endFrequency);
#endif

 private:
  std::unique_ptr<TileBeam> _tileBeam;
  casacore::MPosition _ant1Pos;
  casacore::MEpoch _time;
  double _frequency;
};

#endif
