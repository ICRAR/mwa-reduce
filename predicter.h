#ifndef PREDICTER_H
#define PREDICTER_H

#include <complex>
#include <vector>

#include "aocommon/threadpool.h"

class Predicter {
 public:
  typedef double NumType;
  typedef std::complex<NumType> CNumType;

  Predicter(NumType phaseCentreRA, NumType phaseCentreDec,
            NumType startFrequency, NumType endFrequency, size_t channelCount)
      : _ra0(phaseCentreRA),
        _dec0(phaseCentreDec),
        _startFrequency(startFrequency),
        _endFrequency(endFrequency),
        _channelCount(channelCount) {
    _totalFlux[0] = 0.0;
    _totalFlux[1] = 0.0;
    _totalFlux[2] = 0.0;
    _totalFlux[3] = 0.0;
  }

  /**
   * Initializes the l and m position(s) of the source.
   */
  void Initialize(class ModelSource &source,
                  class BeamEvaluator *beamEvaluator = nullptr);
  void Initialize(class Model &model,
                  class BeamEvaluator *beamEvaluator = nullptr);
  void ReportSources(class Model &model);
  void UpdateBeam(class Model &model, size_t startChannel, size_t endChannel);

  // CNumType Predict(const class ModelSource &source, NumType u, NumType v,
  // NumType w, size_t channelIndex, size_t polarizationIndex); CNumType
  // Predict(const class Model &model, NumType u, NumType v, NumType w, size_t
  // channelIndex, size_t polarizationIndex);

  void Predict4(CNumType *dest, const class ModelSource &source, NumType u,
                NumType v, NumType w, size_t channelIndex, size_t a1,
                size_t a2);
  void Predict4(CNumType *dest, const class Model &model, NumType u, NumType v,
                NumType w, size_t channelIndex, size_t a1, size_t a2);

  NumType TotalFlux(size_t p) { return std::fabs(_totalFlux[p]); }

 private:
#ifndef CUDA_SUPPORT
  void initialize(class ModelComponent &source);
  void updateBeam(class ModelComponent &source, size_t startChannel,
                  size_t endChannel);
#else
  void initialize_one(class ModelComponent &source);
  void initialize_three(std::vector<ModelComponent *> &all_components);
  void updateBeam(std::vector<ModelComponent *> &all_components,
                  size_t startChannel, size_t endChannel);
#endif
  void predict4(CNumType *dest, const class ModelComponent &component,
                NumType u, NumType v, NumType w, size_t channelIndex, size_t a1,
                size_t a2);
  struct SourceParameters {
    NumType l, m, lmsqrt;
    CNumType *brightness, *beamValues, *appBrightness;
    NumType gausTransf[4];
  };

  NumType _ra0, _dec0, _startFrequency, _endFrequency;
  size_t _channelCount;
  class BeamEvaluator *_beamEvaluator;
  CNumType _totalFlux[4];
  ThreadPool _threads;
};

#endif
