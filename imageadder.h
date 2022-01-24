#ifndef IMAGE_ADDER_H
#define IMAGE_ADDER_H

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/uvector.h>
#include "image.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <memory>

class ImageAdder {
 public:
  ImageAdder()
      : _width(0),
        _height(0),
        _count(0),
        _outImage(0),
        _outWeights(0),
        _imgWriter(),
        _frequencySum(0.0),
        _lowestFreq(0.0),
        _highestFreq(0.0),
        _normalizationFactorSum(0.0),
        _normalizationFactorWeightSum(0.0),
        _totalImageWeight(0.0),
        _totalImageWeightWeightSum(0.0),
        _threshold(0.0),
        _weighting(true),
        _square(false) {}

  void Add(const std::string& imageFilename,
           const std::string& weightFilename = "-1",
           const std::string& imaginaryWeightFilename = "") {
    aocommon::FitsReader inpReader(imageFilename);
    if (_outImage.empty()) {
      _width = inpReader.ImageWidth(), _height = inpReader.ImageHeight();
      _imgWriter.reset(new aocommon::FitsWriter(inpReader));

      const size_t size = _width * _height;
      _outImage.assign(size, 0.0);
      _outWeights.assign(size, 0.0);
    } else {
      if (_width != inpReader.ImageWidth() ||
          _height != inpReader.ImageHeight())
        throw std::runtime_error("Not all images have same size");
    }

    double thisFrequency = inpReader.Frequency();
    _frequencySum += thisFrequency;
    if (_count == 0) {
      _lowestFreq = thisFrequency - inpReader.Bandwidth() * 0.5;
      _highestFreq = thisFrequency + inpReader.Bandwidth() * 0.5;
    } else {
      if (thisFrequency - inpReader.Bandwidth() * 0.5 < _lowestFreq)
        _lowestFreq = thisFrequency - inpReader.Bandwidth() * 0.5;
      if (thisFrequency + inpReader.Bandwidth() * 0.5 > _highestFreq)
        _highestFreq = thisFrequency + inpReader.Bandwidth() * 0.5;
    }
    _count++;

    aocommon::UVector<double> inpImage(_width * _height),
        beamImage(_width * _height);

    inpReader.Read<double>(&inpImage[0]);
    double wscImageWeight = 1.0;
    if (_square) {
      for (double& v : inpImage) v = v * v;
    }
    if (_weighting) {
      if (!inpReader.ReadDoubleKeyIfExists("WSCIMGWG", wscImageWeight)) {
        std::cerr << "Warning: Keyword WSCIMGWG not found!\n";
      }
    }
    double wscNormalizationFactor;
    if (!inpReader.ReadDoubleKeyIfExists("WSCNORMF", wscNormalizationFactor)) {
      std::cerr << "Warning: Keyword WSCNORM not found!\n";
      wscNormalizationFactor = 1.0;
    }

    if (std::string(weightFilename) == "-1") {
      beamImage.assign(_width * _height, 1.0);
    } else if (!imaginaryWeightFilename.empty()) {
      aocommon::FitsReader realReader(weightFilename),
          imagReader(imaginaryWeightFilename);
      if (realReader.ImageWidth() != _width ||
          realReader.ImageHeight() != _height)
        throw std::runtime_error("Real beam and image do not have same size");
      if (imagReader.ImageWidth() != _width ||
          imagReader.ImageHeight() != _height)
        throw std::runtime_error(
            "Imaginary beam and image do not have same size");

      aocommon::UVector<double> realImage(_width * _height),
          imagImage(_width * _height);
      realReader.Read<double>(&realImage[0]);
      imagReader.Read<double>(&imagImage[0]);
      for (size_t j = 0; j != _width * _height; ++j) {
        double r = realImage[j], i = imagImage[j];
        beamImage[j] = (r * r + i * i);
      }
    } else {
      aocommon::FitsReader weightsReader(weightFilename);
      if (weightsReader.ImageWidth() != _width ||
          weightsReader.ImageHeight() != _height)
        throw std::runtime_error("Weights and image do not have same size");
      weightsReader.Read<double>(&beamImage[0]);
    }

    bool acceptImage = true;
    if (_threshold != 0.0) {
      double rms = DImage::RMS(inpImage.data(), inpImage.size());
      if (rms == 0.0 || !std::isfinite(rms) || rms > _threshold)
        acceptImage = false;
    }

    if (acceptImage) {
      double centralWeight = beamImage[_width / 2 + (_height / 2) * _width];
      _normalizationFactorSum +=
          wscNormalizationFactor * wscImageWeight * centralWeight;
      _normalizationFactorWeightSum += wscImageWeight * centralWeight;
      _totalImageWeight += wscImageWeight * centralWeight;
      _totalImageWeightWeightSum += centralWeight;

      // Add the images in
      double *outImagePtr = _outImage.data(),
             *outWeightPtr = _outWeights.data();
      aocommon::UVector<double>::iterator inpBeamIter = beamImage.begin();
      for (aocommon::UVector<double>::iterator i = inpImage.begin();
           i != inpImage.end(); ++i) {
        double beamVal = *inpBeamIter;
        *outImagePtr += (*i) * beamVal * wscImageWeight;
        *outWeightPtr += beamVal * beamVal * wscImageWeight;

        ++inpBeamIter;
        ++outImagePtr;
        ++outWeightPtr;
      }
      std::cout << '.' << std::flush;
    } else {
      std::cout << 'R' << std::flush;
    }
  }

  void Finish(const std::string& outputFilename,
              const std::string& outputWeightsFilename = "") {
    // Divide the weight out
    const double* weightsIter = _outWeights.data();
    for (double& value : _outImage) {
      double weight = *weightsIter;
      value /= weight;

      ++weightsIter;
    }
    if (_square) {
      for (double& value : _outImage) {
        value = sqrt(value);
      }
    }

    _imgWriter->SetFrequency(_frequencySum / _count,
                             (_highestFreq - _lowestFreq));
    _imgWriter->SetExtraKeyword(
        "WSCNORMF", _normalizationFactorSum / _normalizationFactorWeightSum);
    _imgWriter->SetExtraKeyword("WSCIMGWG",
                                _totalImageWeight / _totalImageWeightWeightSum);
    _imgWriter->Write<double>(outputFilename, _outImage.data());

    if (!outputWeightsFilename.empty() && outputWeightsFilename != "-")
      _imgWriter->Write<double>(outputWeightsFilename, _outWeights.data());
  }

  void SetThreshold(double threshold) { _threshold = threshold; }

  void SetWeighting(bool weighting) { _weighting = weighting; }

  void SetSquare(bool square) { _square = square; }

  void SetSICorrection(std::vector<double>&& siCorrection) {
    _siCorrection = std::move(siCorrection);
  }

 private:
  size_t _width, _height, _count;
  aocommon::UVector<double> _outImage, _outWeights;
  std::unique_ptr<aocommon::FitsWriter> _imgWriter;
  double _frequencySum, _lowestFreq, _highestFreq;
  double _normalizationFactorSum, _normalizationFactorWeightSum;
  double _totalImageWeight, _totalImageWeightWeightSum;
  double _threshold;
  bool _weighting, _square;
  std::vector<double> _siCorrection;
};

#endif
