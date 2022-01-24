#ifndef IMAGE_H
#define IMAGE_H

#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/uvector.h>

#include <cstring>
#include <cmath>
#include <complex>
#include <memory>

template <typename NumT>
class ImageT {
 public:
  typedef NumT value_type;
  typedef value_type num_t;

  typedef value_type* iterator;
  typedef const value_type* const_iterator;

  typedef std::unique_ptr<ImageT<NumT>> Ptr;

  ImageT() noexcept : _data(nullptr), _width(0), _height(0) {}
  ImageT(size_t width, size_t height);
  ImageT(size_t width, size_t height, value_type initialValue);

  ~ImageT() noexcept;

  ImageT(const ImageT<NumT>&);
  ImageT<NumT>& operator=(const ImageT<NumT>&);

  ImageT<NumT>& operator=(value_type value) noexcept {
    std::fill_n(_data, _width * _height, value);
    return *this;
  }

  ImageT<NumT>& Assign(const NumT* begin, const NumT* end) {
    if (size_t(end - begin) != _width * _height)
      throw std::runtime_error("Invalid Assign()");
    std::copy(begin, end, _data);
    return *this;
  }

  ImageT(ImageT<NumT>&& source) noexcept;
  ImageT<NumT>& operator=(ImageT<NumT>&& source) noexcept;

  static Ptr Make(size_t width, size_t height) {
    return Ptr(new ImageT<NumT>(width, height));
  }
  static Ptr Make(size_t width, size_t height, value_type initialValue) {
    return Ptr(new ImageT<NumT>(width, height, initialValue));
  }

  value_type* data() { return _data; }
  const value_type* data() const { return _data; }

  size_t Width() const { return _width; }
  size_t Height() const { return _height; }
  size_t size() const { return _width * _height; }
  bool empty() const { return _width == 0 || _height == 0; }

  iterator begin() { return _data; }
  const_iterator begin() const { return _data; }

  iterator end() { return _data + _width * _height; }
  const_iterator end() const { return _data + _width * _height; }

  const value_type& operator[](size_t index) const { return _data[index]; }
  value_type& operator[](size_t index) { return _data[index]; }

  ImageT<NumT>& operator+=(const ImageT<NumT>& other);
  ImageT<NumT>& operator-=(const ImageT<NumT>& other);

  ImageT<NumT>& operator*=(value_type factor);
  ImageT<NumT>& operator*=(const ImageT<NumT>& other);
  ImageT<NumT>& operator/=(value_type factor) {
    return (*this) *= value_type(1.0) / factor;
  }

  ImageT<NumT>& Sqrt();
  ImageT<NumT>& Square();
  ImageT<NumT>& SquareWithFactor(NumT factor);
  ImageT<NumT>& AddWithFactor(const ImageT<NumT>& rhs, NumT factor);
  ImageT<NumT>& AddSquared(const ImageT<NumT>& rhs);
  ImageT<NumT>& AddSquared(const ImageT<NumT>& rhs, NumT factor);

  void reset();

  /** Cut-off the borders of an image.
   * @param outWidth Should be &lt;= inWidth.
   * @param outHeight Should be &lt;= inHeight.
   */
  static void Trim(value_type* output, size_t outWidth, size_t outHeight,
                   const value_type* input, size_t inWidth, size_t inHeight);

  ImageT<NumT> Trim(size_t outWidth, size_t outHeight) const {
    ImageT<NumT> image(outWidth, outHeight);
    Trim(image.data(), outWidth, outHeight, data(), Width(), Height());
    return image;
  }

  ImageT<NumT> TrimBox(size_t x1, size_t y1, size_t boxWidth,
                       size_t boxHeight) const {
    ImageT<NumT> image(boxWidth, boxHeight);
    TrimBox(image.data(), x1, y1, boxWidth, boxHeight, data(), Width(),
            Height());
    return image;
  }

  template <typename T>
  static void TrimBox(T* output, size_t x1, size_t y1, size_t boxWidth,
                      size_t boxHeight, const T* input, size_t inWidth,
                      size_t inHeight);

  template <typename T>
  static void CopyMasked(T* to, size_t toX, size_t toY, size_t toWidth,
                         const T* from, size_t fromWidth, size_t fromHeight,
                         const bool* fromMask);

  /** Extend an image with zeros, complement of Trim.
   * @param outWidth Should be &gt;= inWidth.
   * @param outHeight Should be &gt;= inHeight.
   */
  static void Untrim(value_type* output, size_t outWidth, size_t outHeight,
                     const value_type* input, size_t inWidth, size_t inHeight);

  ImageT<NumT> Untrim(size_t outWidth, size_t outHeight) const {
    ImageT<NumT> image(outWidth, outHeight);
    Untrim(image.data(), outWidth, outHeight, data(), Width(), Height());
    return image;
  }

  static value_type Median(const value_type* data, size_t size) {
    aocommon::UVector<value_type> copy;
    return median_with_copy(data, size, copy);
  }

  static value_type MAD(const value_type* data, size_t size);

  value_type Sum() const;
  value_type Average() const;

  value_type Min() const;
  value_type Max() const;

  value_type StdDevFromMAD() const {
    return StdDevFromMAD(_data, _width * _height);
  }
  static value_type StdDevFromMAD(const value_type* data, size_t size) {
    // norminv(0.75) x MAD
    return value_type(1.48260221850560) * MAD(data, size);
  }

  value_type RMS() const { return RMS(_data, _width * _height); }
  static value_type RMS(const value_type* data, size_t size) {
    value_type sum = 0.0;
    for (size_t i = 0; i != size; ++i) sum += data[i] * data[i];
    return std::sqrt(sum / value_type(size));
  }

  void Negate() {
    for (value_type& d : *this) d = -d;
  }

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);

 private:
  value_type* _data;
  size_t _width, _height;

  static value_type median_with_copy(const value_type* data, size_t size,
                                     aocommon::UVector<value_type>& copy);
};

template <typename NumT>
template <typename T>
void ImageT<NumT>::TrimBox(T* output, size_t x1, size_t y1, size_t boxWidth,
                           size_t boxHeight, const T* input, size_t inWidth,
                           size_t /*inHeight*/) {
  size_t endY = y1 + boxHeight;
  for (size_t y = y1; y != endY; ++y) {
    std::copy_n(&input[y * inWidth + x1], boxWidth,
                &output[(y - y1) * boxWidth]);
  }
}

template <typename NumT>
template <typename T>
void ImageT<NumT>::CopyMasked(T* to, size_t toX, size_t toY, size_t toWidth,
                              const T* from, size_t fromWidth,
                              size_t fromHeight, const bool* fromMask) {
  for (size_t y = 0; y != fromHeight; ++y) {
    for (size_t x = 0; x != fromWidth; ++x) {
      if (fromMask[y * fromWidth + x])
        to[toX + (toY + y) * toWidth + x] = from[y * fromWidth + x];
    }
  }
}

typedef ImageT<double> DImage;
typedef ImageT<float> Image;
typedef ImageT<std::complex<float>> ImageCF;

template <class T>
using ImageTC = ImageT<std::complex<T>>;

#endif
