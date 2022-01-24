
#include "modelrenderer.h"

#include "../model/model.h"

#include "../fftconvolver.h"
#include "../fftwmanager.h"

#include <aocommon/imagecoordinates.h>
#include <aocommon/uvector.h>

#include <cmath>

#include <boost/algorithm/clamp.hpp>

using aocommon::ImageCoordinates;
using boost::algorithm::clamp;

template <typename T>
T ModelRenderer::gaus(T x, T sigma) {
  long double xi = x / sigma;
  return exp(T(-0.5) * xi * xi);  // / (sigma * sqrt(T(2.0) * M_PIl));
}

/** Restore a circular beam*/
void ModelRenderer::Restore(double* imageData, size_t imageWidth,
                            size_t imageHeight, const Model& model,
                            long double beamSize, long double startFrequency,
                            long double endFrequency,
                            aocommon::PolarizationEnum polarization) {
  // Using the FWHM formula for a Gaussian:
  long double sigma = beamSize / (2.0L * sqrtl(2.0L * logl(2.0L)));

  int boundingBoxSize =
      ceil(sigma * 20.0 / std::min(_pixelScaleL, _pixelScaleM));
  for (const ModelSource& src : model) {
    for (const ModelComponent& comp : src) {
      const long double posRA = comp.PosRA();
      const long double posDec = comp.PosDec();
      long double sourceL;
      long double sourceM;
      ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA,
                                  _phaseCentreDec, sourceL, sourceM);
      const SpectralEnergyDistribution& sed = comp.SED();
      const long double intFlux =
          sed.IntegratedFlux(startFrequency, endFrequency, polarization);

      int sourceX, sourceY;
      ImageCoordinates::LMToXY<long double>(
          sourceL - _phaseCentreDL, sourceM - _phaseCentreDM, _pixelScaleL,
          _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);

      const int xLeft = clamp(sourceX - boundingBoxSize, 0, int(imageWidth));
      const int xRight =
          clamp(sourceX + boundingBoxSize, xLeft, int(imageWidth));
      const int yTop = clamp(sourceY - boundingBoxSize, 0, int(imageHeight));
      const int yBottom =
          clamp(sourceY + boundingBoxSize, yTop, int(imageHeight));

      for (int y = yTop; y != yBottom; ++y) {
        double* imageDataPtr = imageData + y * imageWidth + xLeft;
        for (int x = xLeft; x != xRight; ++x) {
          long double l, m;
          ImageCoordinates::XYToLM<long double>(
              x, y, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, l, m);
          l += _phaseCentreDL;
          m += _phaseCentreDM;
          long double dist = sqrt((l - sourceL) * (l - sourceL) +
                                  (m - sourceM) * (m - sourceM));
          long double g = gaus(dist, sigma);
          (*imageDataPtr) += double(g * intFlux);
          ++imageDataPtr;
        }
      }
    }
  }
}

void ModelRenderer::RenderGaussianComponent(
    float* imageData, size_t imageWidth, size_t imageHeight,
    long double phaseCentreRA, long double phaseCentreDec,
    long double pixelScaleL, long double pixelScaleM, long double phaseCentreDL,
    long double phaseCentreDM, long double posRA, long double posDec,
    long double gausMaj, long double gausMin, long double gausPA,
    long double flux) {
  // Using the FWHM formula for a Gaussian:
  const long double fwhm_constant = (2.0L * std::sqrt(2.0L * std::log(2.0L)));
  const long double sigmaMaj = gausMaj / fwhm_constant;
  const long double sigmaMin = gausMin / fwhm_constant;
  // TODO this won't work for non-equally spaced dimensions
  const long double minPixelScale = std::min(pixelScaleL, pixelScaleM);

  // Make rotation matrix
  long double transf[4];
  // Position angle is angle from North:
  const long double angle = gausPA + 0.5 * M_PI;
  transf[2] = std::sin(angle);
  transf[0] = std::cos(angle);
  transf[1] = -transf[2];
  transf[3] = transf[0];
  const double sigmaMax = std::max(std::fabs(sigmaMaj * transf[0]),
                                   std::fabs(sigmaMaj * transf[1]));
  // Multiply with scaling matrix to make variance 1.
  transf[0] = transf[0] / sigmaMaj;
  transf[1] = transf[1] / sigmaMaj;
  transf[2] = transf[2] / sigmaMin;
  transf[3] = transf[3] / sigmaMin;
  const int boundingBoxSize = std::ceil(sigmaMax * 20.0 / minPixelScale);
  long double sourceL, sourceM;
  ImageCoordinates::RaDecToLM(posRA, posDec, phaseCentreRA, phaseCentreDec,
                              sourceL, sourceM);

  // Calculate the bounding box
  int sourceX, sourceY;
  ImageCoordinates::LMToXY<long double>(
      sourceL - phaseCentreDL, sourceM - phaseCentreDM, pixelScaleL,
      pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
  const int xLeft = clamp(sourceX - boundingBoxSize, 0, int(imageWidth));
  const int xRight = clamp(sourceX + boundingBoxSize, xLeft, int(imageWidth));
  const int yTop = clamp(sourceY - boundingBoxSize, 0, int(imageHeight));
  const int yBottom = clamp(sourceY + boundingBoxSize, yTop, int(imageHeight));

  aocommon::UVector<double> values;
  double fluxSum = 0.0;
  for (int y = yTop; y != yBottom; ++y) {
    for (int x = xLeft; x != xRight; ++x) {
      long double l, m;
      ImageCoordinates::XYToLM<long double>(x, y, pixelScaleL, pixelScaleM,
                                            imageWidth, imageHeight, l, m);
      l += phaseCentreDL;
      m += phaseCentreDM;
      const long double lTransf =
          (l - sourceL) * transf[0] + (m - sourceM) * transf[1];
      const long double mTransf =
          (l - sourceL) * transf[2] + (m - sourceM) * transf[3];
      const long double dist = std::sqrt(lTransf * lTransf + mTransf * mTransf);
      long double v = gaus(dist, 1.0L) * flux;
      fluxSum += double(v);
      values.emplace_back(v);
    }
  }
  const double* iter = values.data();
  const double factor = flux / fluxSum;
  for (int y = yTop; y != yBottom; ++y) {
    float* imageDataPtr = imageData + y * imageWidth + xLeft;
    for (int x = xLeft; x != xRight; ++x) {
      (*imageDataPtr) += *iter * factor;
      ++imageDataPtr;
      ++iter;
    }
  }
}

void ModelRenderer::renderPointComponent(float* imageData, size_t imageWidth,
                                         size_t imageHeight, long double posRA,
                                         long double posDec, long double flux) {
  long double sourceL, sourceM;
  ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec,
                              sourceL, sourceM);
  sourceL -= _phaseCentreDL;
  sourceM -= _phaseCentreDM;

  int sourceX, sourceY;
  ImageCoordinates::LMToXY<long double>(sourceL, sourceM, _pixelScaleL,
                                        _pixelScaleM, imageWidth, imageHeight,
                                        sourceX, sourceY);

  if (sourceX >= 0 && sourceX < (int)imageWidth && sourceY >= 0 &&
      sourceY < (int)imageHeight) {
    float* imageDataPtr = imageData + sourceY * imageWidth + sourceX;
    (*imageDataPtr) += double(flux);
  }
}

/** Restore a model with an elliptical beam */
void ModelRenderer::Restore(float* imageData, size_t imageWidth,
                            size_t imageHeight, const Model& model,
                            long double beamMaj, long double beamMin,
                            long double beamPA, long double startFrequency,
                            long double endFrequency,
                            aocommon::PolarizationEnum polarization,
                            size_t threadCount) {
  aocommon::UVector<float> renderedWithoutBeam(imageWidth * imageHeight, 0.0);
  RenderModel(renderedWithoutBeam.data(), imageWidth, imageHeight, model,
              startFrequency, endFrequency, polarization);
  Restore(imageData, renderedWithoutBeam.data(), imageWidth, imageHeight,
          beamMaj, beamMin, beamPA, _pixelScaleL, _pixelScaleM, threadCount);
}

/**
 * Restore a diffuse image (e.g. produced with multi-scale clean)
 */
void ModelRenderer::Restore(float* imageData, const float* modelData,
                            size_t imageWidth, size_t imageHeight,
                            long double beamMaj, long double beamMin,
                            long double beamPA, long double pixelScaleL,
                            long double pixelScaleM, size_t threadCount) {
  if (beamMaj == 0.0 && beamMin == 0.0) {
    for (size_t j = 0; j != imageWidth * imageHeight; ++j) {
      imageData[j] += modelData[j];
    }
  } else {
    // Using the FWHM formula for a Gaussian:
    long double sigmaMaj = beamMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
    long double sigmaMin = beamMin / (2.0L * sqrtl(2.0L * logl(2.0L)));

    // Make rotation matrix
    long double transf[4];
    // Position angle is angle from North:
    long double angle = beamPA + 0.5 * M_PI;
    transf[2] = sin(angle);
    transf[0] = cos(angle);
    transf[1] = -transf[2];
    transf[3] = transf[0];
    double sigmaMax = std::max(std::fabs(sigmaMaj * transf[0]),
                               std::fabs(sigmaMaj * transf[1]));
    // Multiple with scaling matrix to make variance 1.
    transf[0] = transf[0] / sigmaMaj;
    transf[1] = transf[1] / sigmaMaj;
    transf[2] = transf[2] / sigmaMin;
    transf[3] = transf[3] / sigmaMin;

    size_t minDimension = std::min(imageWidth, imageHeight);
    size_t boundingBoxSize = std::min<size_t>(
        ceil(sigmaMax * 40.0 / std::min(pixelScaleL, pixelScaleM)),
        minDimension);
    if (boundingBoxSize % 2 != 0) {
      ++boundingBoxSize;
    }
    if (boundingBoxSize > std::min(imageWidth, imageHeight))
      boundingBoxSize = std::min(imageWidth, imageHeight);
    aocommon::UVector<float> kernel(boundingBoxSize * boundingBoxSize);
    auto iter = kernel.begin();
    for (size_t y = 0; y != boundingBoxSize; ++y) {
      for (size_t x = 0; x != boundingBoxSize; ++x) {
        long double l, m;
        ImageCoordinates::XYToLM<long double>(x, y, pixelScaleL, pixelScaleM,
                                              boundingBoxSize, boundingBoxSize,
                                              l, m);
        long double lTransf = l * transf[0] + m * transf[1],
                    mTransf = l * transf[2] + m * transf[3];
        long double dist = std::sqrt(lTransf * lTransf + mTransf * mTransf);
        *iter = gaus(dist, (long double)1.0);
        ++iter;
      }
    }

    aocommon::UVector<float> convolvedModel(
        modelData, modelData + imageWidth * imageHeight);

    FFTWManager fftw;
    FFTConvolver::Convolve(fftw, convolvedModel.data(), imageWidth, imageHeight,
                           kernel.data(), boundingBoxSize, threadCount);
    for (size_t j = 0; j != imageWidth * imageHeight; ++j) {
      imageData[j] += convolvedModel[j];
    }
  }
}

/**
 * Render without beam convolution, such that each point-source is one pixel.
 */
void ModelRenderer::RenderModel(float* imageData, size_t imageWidth,
                                size_t imageHeight, const Model& model,
                                long double startFrequency,
                                long double endFrequency,
                                aocommon::PolarizationEnum polarization) {
  for (Model::const_iterator src = model.begin(); src != model.end(); ++src) {
    for (ModelSource::const_iterator comp = src->begin(); comp != src->end();
         ++comp) {
      const long double posRA = comp->PosRA(), posDec = comp->PosDec(),
                        intFlux = comp->SED().IntegratedFlux(
                            startFrequency, endFrequency, polarization);

      if (comp->Type() == ModelComponent::GaussianSource) {
        long double gausMaj = comp->MajorAxis(), gausMin = comp->MinorAxis(),
                    gausPA = comp->PositionAngle();
        RenderGaussianComponent(
            imageData, imageWidth, imageHeight, _phaseCentreRA, _phaseCentreDec,
            _pixelScaleL, _pixelScaleM, _phaseCentreDL, _phaseCentreDM, posRA,
            posDec, gausMaj, gausMin, gausPA, intFlux);
      } else
        renderPointComponent(imageData, imageWidth, imageHeight, posRA, posDec,
                             intFlux);
    }
  }
}
