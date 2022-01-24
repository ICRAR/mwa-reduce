#ifndef INTERPOLATED_RENDERER_H
#define INTERPOLATED_RENDERER_H

#include <cstring>

#include <aocommon/uvector.h>

class InterpolatingRenderer {
 public:
  InterpolatingRenderer(size_t kernelSize)
      : _kernelSize(kernelSize + (kernelSize + 1) % 2),
        _hKernel(_kernelSize),
        _vKernel(_kernelSize) {}

  /**
   * This will render a source and sinc-interpolate it so it
   * can be on non-integer positions.
   */
  void RenderSource(float* image, size_t width, size_t height, float flux,
                    double x, double y);

  /**
   * This will render a source and sinc-interpolate with a window it so it
   * can be on non-integer positions.
   */
  void RenderWindowedSource(float* image, size_t width, size_t height,
                            float flux, float x, float y);

 private:
  size_t _kernelSize;
  aocommon::UVector<float> _hKernel;
  aocommon::UVector<float> _vKernel;
};

#endif
