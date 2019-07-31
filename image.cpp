#include "image.h"

#include <algorithm>
#include <cmath>

Image::Image(size_t width, size_t height, ImageBufferAllocator& allocator) :
	_data(allocator.Allocate(width*height)),
	_width(width), _height(height), _allocator(&allocator)
{
}

Image::Image(size_t width, size_t height, double initialValue, ImageBufferAllocator& allocator) :
	_data(allocator.Allocate(width*height)),
	_width(width), _height(height), _allocator(&allocator)
{
	std::fill(_data, _data+_width*_height, initialValue);
}

Image::~Image()
{
	if(_allocator != nullptr)
		_allocator->Free(_data);
}

Image::Image(const Image& source) :
	_data((source._allocator==nullptr) ? nullptr : source._allocator->Allocate(source._width * source._height)),
	_width(source._width),
	_height(source._height),
	_allocator(source._allocator)
{
	std::copy(source._data, source._data + _width*_height, _data);
}

Image& Image::operator=(const Image& source)
{
	if(source._allocator == nullptr)
	{
		reset();
	}
	else {
		if(_allocator == nullptr)
		{
			_allocator = source._allocator;
			_data = _allocator->Allocate(source._width * source._height);
		}
		else if(_width * _height != source._width * source._height || source._allocator != _allocator)
		{
			_allocator->Free(_data);
			_allocator = source._allocator;
			_data = _allocator->Allocate(source._width * source._height);
		}
		_width = source._width;
		_height = source._height;
		std::copy(source._data, source._data + _width*_height, _data);
	}
	return *this;
}

Image& Image::operator=(double value)
{
	for(double& v : *this)
		v = value;
	return *this;
}

Image::Image(Image&& source) :
	_data(source._data),
	_width(source._width),
	_height(source._height),
	_allocator(source._allocator)
{
	source._width = 0;
	source._height = 0;
	source._data = nullptr;
	source._allocator = nullptr;
}

Image& Image::operator=(Image&& source)
{
	std::swap(_data, source._data);
	std::swap(_width, source._width);
	std::swap(_height, source._height);
	std::swap(_allocator, source._allocator);
	return *this;
}

void Image::reset()
{
	if(_allocator != nullptr)
		_allocator->Free(_data);
	_data = nullptr;
	_width = 0;
	_height = 0;
	_allocator = nullptr;
}

Image& Image::operator+=(const Image& other)
{
	for(size_t i=0; i!=_width*_height; ++i)
		_data[i] += other[i];
	return *this;
}

Image& Image::operator*=(double factor)
{
	for(size_t i=0; i!=_width*_height; ++i)
		_data[i] *= factor;
	return *this;
}

Image& Image::operator*=(const Image& other)
{
	for(size_t i=0; i!=_width*_height; ++i)
		_data[i] *= other[i];
	return *this;
}

// Cut-off the borders of an image.
// @param outWidth Should be <= inWidth.
// @param outHeight Should be <= inHeight.
void Image::Trim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight)
{
	size_t startX = (inWidth - outWidth) / 2;
	size_t startY = (inHeight - outHeight) / 2;
	size_t endY = (inHeight + outHeight) / 2;
	for(size_t y=startY; y!=endY; ++y)
	{
		memcpy(&output[(y-startY)*outWidth], &input[y*inWidth + startX], outWidth*sizeof(double));
	}
}

void Image::TrimBox(bool* output, size_t x1, size_t y1, size_t boxWidth, size_t boxHeight, const bool* input, size_t inWidth, size_t /*inHeight*/)
{
	size_t endY = y1 + boxHeight;
	for(size_t y=y1; y!=endY; ++y)
	{
		memcpy(&output[(y-y1)*boxWidth], &input[y*inWidth + x1], boxWidth*sizeof(bool));
	}
}

/** Extend an image with zeros, complement of Trim.
	* @param outWidth Should be &gt;= inWidth.
	* @param outHeight Should be &gt;= inHeight.
	*/
void Image::Untrim(double* output, size_t outWidth, size_t outHeight, const double* input, size_t inWidth, size_t inHeight)
{
	size_t startX = (outWidth - inWidth) / 2;
	size_t endX = (outWidth + inWidth) / 2;
	size_t startY = (outHeight - inHeight) / 2;
	size_t endY = (outHeight + inHeight) / 2;
	for(size_t y=0; y!=startY; ++y)
	{
		double* ptr = &output[y*outWidth];
		for(size_t x=0; x!=outWidth; ++x)
			ptr[x] = 0.0;
	}
	for(size_t y=startY; y!=endY; ++y)
	{
		double* ptr = &output[y*outWidth];
		for(size_t x=0; x!=startX; ++x)
			ptr[x] = 0.0;
		memcpy(&output[y*outWidth + startX], &input[(y-startY)*inWidth], inWidth*sizeof(double));
		for(size_t x=endX; x!=outWidth; ++x)
			ptr[x] = 0.0;
	}
	for(size_t y=endY; y!=outHeight; ++y)
	{
		double* ptr = &output[y*outWidth];
		for(size_t x=0; x!=outWidth; ++x)
			ptr[x] = 0.0;
	}
}

double Image::Sum() const
{
	double sum = 0.0;
	for(const double& v : *this)
		sum += v;
	return sum;
}

double Image::Average() const
{
	return Sum() / size();
}

double Image::Min() const
{
	return *std::min_element(begin(), end());
}

double Image::Max() const
{
	return *std::max_element(begin(), end());
}

double Image::median_with_copy(const double* data, size_t size, ao::uvector<double>& copy)
{
	copy.reserve(size);
	for(const double* i=data ; i!=data+size; ++i)
	{
		if(std::isfinite(*i))
			copy.push_back(*i);
	}
	if(copy.empty())
		return 0.0;
	else {
		bool even = (copy.size()%2) == 0;
		ao::uvector<double>::iterator mid = copy.begin()+(copy.size()-1)/2;
		std::nth_element(copy.begin(), mid, copy.end());
		double median = *mid;
		if(even)
		{
			std::nth_element(mid, mid+1, copy.end());
			median = (median + *(mid+1)) * 0.5;
		}
		return median;
	}
}

double Image::MAD(const double* data, size_t size)
{
	ao::uvector<double> copy;
	double median = median_with_copy(data, size, copy);
	if(copy.empty())
		return 0.0;
		
	// Replace all values by the difference from the mean
	ao::uvector<double>::iterator mid = copy.begin()+(copy.size()-1)/2;
	for(ao::uvector<double>::iterator i=copy.begin(); i!=mid+1; ++i)
		*i = median - *i;
	for(ao::uvector<double>::iterator i=mid+1; i!=copy.end(); ++i)
		*i = *i - median;
	
	std::nth_element(copy.begin(), mid, copy.end());
	median = *mid;
	bool even = (copy.size()%2) == 0;
	if(even)
	{
		std::nth_element(mid, mid+1, copy.end());
		median = (median + *(mid+1)) * 0.5;
	}
	return median;
}
