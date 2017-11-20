#ifndef IMAGE_H
#define IMAGE_H

#include <cstring>

#include "wsclean/imagebufferallocator.h"

#include "uvector.h"

class Image
{
public:
	typedef double* iterator;
	typedef const double* const_iterator;
	
	Image() : _data(nullptr), _width(0), _height(0), _allocator(nullptr) { }
	Image(size_t width, size_t height, class ImageBufferAllocator& allocator);
	Image(size_t width, size_t height, double initialValue, ImageBufferAllocator& allocator);
	
	~Image();
	
	Image(const Image&);
	Image& operator=(const Image&);
	Image& operator=(double value);
	
	Image(Image&& source);
	Image& operator=(Image&& source);
	
	double* data() { return _data; }
	const double* data() const { return _data; }
	
	size_t Width() const { return _width; }
	size_t Height() const { return _height; }
	size_t size() const { return _width * _height; }
	bool empty() const { return _width == 0 || _height == 0; }
	ImageBufferAllocator& Allocator() const { return *_allocator; }
	
	iterator begin() { return _data; }
	const_iterator begin() const { return _data; }
	
	iterator end() { return _data + _width*_height; }
	const_iterator end() const { return _data + _width*_height; }
	
	const double& operator[](size_t index) const { return _data[index]; }
	double& operator[](size_t index) { return _data[index]; }
	
	Image& operator*=(double factor);
	Image& operator*=(const Image& other);
	
	void reset();
	
	/** Cut-off the borders of an image.
	 * @param outWidth Should be &lt;= inWidth.
	 * @param outHeight Should be &lt;= inHeight.
	 */
	static void Trim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static void TrimBox(bool* output, size_t x1, size_t y1, size_t boxWidth, size_t boxHeight, const bool* input, size_t inWidth, size_t inHeight);
	
	/** Extend an image with zeros, complement of Trim.
	 * @param outWidth Should be &gt;= inWidth.
	 * @param outHeight Should be &gt;= inHeight.
	 */
	static void Untrim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight);
	
	static double Median(const double* data, size_t size)
	{
		ao::uvector<double> copy;
		return median_with_copy(data, size, copy);
	}
	
	static double MAD(const double* data, size_t size);
	
	double Sum() const;
	double Average() const;
	
	double Min() const;
	double Max() const;
	
	double StdDevFromMAD() const { return StdDevFromMAD(_data, _width*_height); }
	static double StdDevFromMAD(const double* data, size_t size)
	{
		// norminv(0.75) x MAD
		return 1.48260221850560 * MAD(data, size);
	}
	
	static double RMS(const double* data, size_t size)
	{
		double sum = 0.0;
		for(size_t i=0; i!=size; ++i)
			sum += data[i]*data[i];
		return sqrt(sum/size);
	}
	
	void Negate()
	{
		for(double& d : *this)
			d = -d;
	}
private:
	double* _data;
	size_t _width, _height;
	ImageBufferAllocator* _allocator;
	
	static double median_with_copy(const double* data, size_t size, ao::uvector<double>& copy);
};

#endif
