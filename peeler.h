#ifndef PEELER_H
#define PEELER_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <mutex>

#include <queue>

#include <aocommon/banddata.h>
#include "lane.h"

#include "model/model.h"

class Peeler {
 public:
  Peeler(casacore::MeasurementSet &ms, size_t threadCount);

  void Perform();

  void SetBeamOnSource(bool beamOnSource) { _beamOnSource = beamOnSource; }

  void SetMWAPath(const std::string &path) { _mwaPath = path; }

  void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }

  void SetOnlyScalar(bool onlyScalar) { _onlyScalar = onlyScalar; }

  void SetOnlyDiag(bool onlyDiag) { _onlyDiag = onlyDiag; }

  void SetOnlyRotation(bool onlyRotation) { _onlyRotation = onlyRotation; }

  void SetModelFilename(const std::string &modelFilename) {
    _modelFilename = modelFilename;
  }

  void SetModel(const Model &model) { _model = model; }

  void SetRHSSolutionFile(const std::string &rhsSolutionFile) {
    _rhsSolutionFile = rhsSolutionFile;
  }

  void SetDataColumnName(const std::string &dataColumnName) {
    _dataColumnName = dataColumnName;
  }

  void SetNIter(size_t nIter) { _nIter = nIter; }

  void SetAccuracy(double minAccuracy, double stoppingAccuracy) {
    _minAccuracy = minAccuracy;
    _stoppingAccuracy = stoppingAccuracy;
  }

  void SetMinUVW(double minUVW) { _minUVW = minUVW; }

  void SetSolutionInterval(size_t solutionInterval) {
    _solutionInterval = solutionInterval;
  }

  void SetSaveSolutionFiles(bool saveSolutionFiles) {
    _saveSolutionsFiles = saveSolutionFiles;
  }

 private:
  struct ThreadData {
    ThreadData() {}

    std::mutex *mutex;
    std::queue<size_t> *tasks;
    std::vector<class CalibrationMethod *> *calMethods;
  };
  struct SubtractThreadInfo {
    size_t rowIndex;
    casacore::Array<casacore::Complex> *data;
    bool readyForWrite;
    double u, v, w;
    size_t a1, a2;
  };

  void calibrateThreadFunction(ThreadData data);

  casacore::MeasurementSet &_ms;
  aocommon::BandData _bandData;

  bool _beamOnSource, _applyBeam, _onlyScalar, _onlyDiag, _onlyRotation,
      _saveSolutionsFiles;
  std::string _modelFilename, _rhsSolutionFile, _dataColumnName, _mwaPath;
  size_t _nIter;
  double _minAccuracy, _stoppingAccuracy, _minUVW;
  size_t _solutionInterval;
  Model _model;
  size_t _threadCount;
};

#endif
