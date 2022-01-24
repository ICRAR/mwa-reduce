#include "interpolatingrenderer.h"

#include <aocommon/uvector.h>

#include <cmath>

void InterpolatingRenderer::RenderSource(float* image, size_t width,
                                         size_t height, float flux, double x,
                                         double y) {
  aocommon::UVector<float> hSinc(std::min<size_t>(width, 128) + 1),
      vSinc(std::min<size_t>(height, 128) + 1);

  int midH = hSinc.size() / 2;
  int midV = vSinc.size() / 2;

  double xr = x - std::floor(x);
  double yr = y - std::floor(y);

  for (size_t i = 0; i != hSinc.size(); ++i) {
    double xi = (int(i) - midH - xr) * M_PI;
    if (xi == 0.0) {
      hSinc[i] = 1.0;
    } else {
      hSinc[i] = std::sin(xi) / xi;
    }
  }
  for (size_t i = 0; i != vSinc.size(); ++i) {
    double yi = (int(i) - midV - yr) * M_PI;
    if (yi == 0.0) {
      vSinc[i] = 1.0;
    } else {
      vSinc[i] = std::sin(yi) / yi;
    }
  }

  const int xOffset = std::floor(x) - midH;
  const int yOffset = std::floor(y) - midV;
  const size_t startX = std::max<int>(xOffset, 0);
  const size_t startY = std::max<int>(yOffset, 0);
  const size_t endX = std::min<size_t>(xOffset + _hKernel.size(), width);
  const size_t endY = std::min<size_t>(yOffset + _vKernel.size(), height);

  for (size_t yi = startY; yi != endY; ++yi) {
    float* ptr = &image[yi * width];
    const double vFlux = flux * vSinc[yi - yOffset];
    for (size_t xi = startX; xi != endX; ++xi) {
      ptr[xi] += vFlux * hSinc[xi - xOffset];
    }
  }
}

void InterpolatingRenderer::RenderWindowedSource(float* image, size_t width,
                                                 size_t height, float flux,
                                                 float x, float y) {
  const size_t n = _hKernel.size();
  const int midH = n / 2;
  const int midV = n / 2;

  const float xr = x - std::floor(x);
  const float yr = y - std::floor(y);

  for (size_t i = 0; i != n; ++i) {
    float xi = (int(i) - midH - xr) * M_PI;
    if (xi == 0.0) {
      _hKernel[i] = 1.0;
    } else {
      const float hannFac = std::cos(xi / n);
      _hKernel[i] = hannFac * hannFac * std::sin(xi) / xi;
    }
  }

  for (size_t i = 0; i != _vKernel.size(); ++i) {
    float yi = (int(i) - midV - yr) * M_PI;
    if (yi == 0.0) {
      _vKernel[i] = 1.0;
    } else {
      const float hannFac = std::cos(yi / n);
      _vKernel[i] = hannFac * hannFac * std::sin(yi) / yi;
    }
  }

  const int xOffset = std::floor(x) - midH;
  const int yOffset = std::floor(y) - midV;
  const size_t startX = std::max<int>(xOffset, 0);
  const size_t startY = std::max<int>(yOffset, 0);
  const size_t endX =
      std::max<int>(std::min<int>(xOffset + int(n), int(width)), int(startX));
  const size_t endY =
      std::max<int>(std::min<int>(yOffset + int(n), int(height)), int(startY));
  for (size_t yi = startY; yi != endY; ++yi) {
    float* ptr = &image[yi * width];
    const float vFlux = flux * _vKernel[yi - yOffset];
    for (size_t xi = startX; xi != endX; ++xi) {
      ptr[xi] += vFlux * _hKernel[xi - xOffset];
    }
  }
}
