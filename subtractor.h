#ifndef SUBTRACTOR_H
#define SUBTRACTOR_H

#include "beamevaluator.h"
#include <aocommon/banddata.h>
#include "predicter.h"
#include "mspredicter.h"
#include "progressbar.h"

#include "model/model.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <cmath>
#include <fstream>

class Subtractor {
 public:
  Subtractor(size_t threadCount)
      : _revert(false),
        _setToModel(false),
        _addNoise(false),
        _applyBeam(false),
        _noiseSigma(1.0),
        _threadCount(threadCount),
        _dataColumn("DATA") {}

  void SetRevert(bool revert) { _revert = revert; }
  void SetToModel(bool setToModel) { _setToModel = setToModel; }
  void SetAddNoise(bool addNoise) { _addNoise = addNoise; }
  void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
  void SetNoiseSigma(double noiseSigma) { _noiseSigma = noiseSigma; }
  void SetDataColumn(const std::string &dataColumn) {
    _dataColumn = dataColumn;
  }
  void SetMWAPath(const std::string &path) { _mwaPath = path; }

  void Subtract(casacore::MeasurementSet &ms, const Model &model) {
    /**
     * Read some meta data from the measurement set
     */
    aocommon::BandData bandData(ms.spectralWindow());
    size_t channelCount = bandData.ChannelCount();

    typedef float num_t;
    typedef std::complex<num_t> complex_t;
    casacore::ArrayColumn<complex_t> dataColumn(ms, _dataColumn);

    casacore::IPosition dataShape = dataColumn.shape(0);
    unsigned polarizationCount = dataShape[0];

    MSPredicter predicter(ms, _threadCount, model);
    predicter.SetApplyBeam(_applyBeam);
    predicter.SetMWAPath(_mwaPath);
    predicter.Start(true);

    if (_addNoise) std::cout << "Adding noise of " << _noiseSigma << " Jy.\n";

    /**
     * Subtract
     */
    std::ostringstream taskDesc;
    if (_revert)
      taskDesc << "Adding back ";
    else if (_setToModel)
      taskDesc << "Setting to model from ";
    else
      taskDesc << "Subtracting ";
    taskDesc << model.SourceCount() << " sources";
    ProgressBar progress(taskDesc.str());

    casacore::Array<complex_t> data(dataShape);
    MSPredicter::RowData rowData;
    while (predicter.GetNextRow(rowData)) {
      size_t rowIndex = rowData.rowIndex;

      std::unique_lock<std::mutex> lock(predicter.IOMutex());
      progress.SetProgress(rowIndex, ms.nrow());
      dataColumn.get(rowIndex, data);
      lock.unlock();

      casacore::Array<complex_t>::iterator dataPtr = data.begin();
      std::complex<double> *modelDataPtr = rowData.modelData;
      for (size_t ch = 0; ch != channelCount; ++ch) {
        for (size_t p = 0; p != polarizationCount; ++p) {
          std::complex<double> predicted;
          if (_revert || _setToModel)
            predicted = *modelDataPtr;
          else
            predicted = -*modelDataPtr;
          if (_addNoise) addGausNoise(predicted, _noiseSigma);
          if (_setToModel)
            *dataPtr = predicted;
          else
            *dataPtr += predicted;
          ++dataPtr;
          ++modelDataPtr;
        }
      }

      lock.lock();
      dataColumn.put(rowIndex, data);
      lock.unlock();

      predicter.FinishRow(rowData);
    }
  }

 private:
  template <typename T>
  void addGausNoise(std::complex<T> &value, double sigma) {
    long double x1, x2, w;

    do {
      long double r1 = (long double)rand() / (long double)RAND_MAX;
      long double r2 = (long double)rand() / (long double)RAND_MAX;
      x1 = 2.0 * r1 - 1.0;
      x2 = 2.0 * r2 - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * std::log(w)) / w) * sigma;
    value += std::complex<T>(x1 * w, x2 * w);
  }

  bool _revert, _setToModel, _addNoise, _applyBeam;
  double _noiseSigma;
  size_t _threadCount;
  std::string _dataColumn, _mwaPath;
};

#endif
