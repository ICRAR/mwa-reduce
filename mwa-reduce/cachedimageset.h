#ifndef CACHED_IMAGE_SET_H
#define CACHED_IMAGE_SET_H

#include "fitswriter.h"
#include "imagebufferallocator.h"
#include "fitsreader.h"

#include <string.h>
#include <set>

class CachedImageSet
{
	typedef double value_t;
public:
	CachedImageSet() : _allocator(0), _image(0)
	{
	}
	
	~CachedImageSet()
	{
		if(_allocator != 0)
			_allocator->Free(_image);
	}
	
	void Initialize(FitsWriter& writer, size_t imageCount, const std::string& prefix, ImageBufferAllocator<value_t>& allocator)
	{
		_writer = writer;
		_imageCount = imageCount;
		_prefix = prefix;
		if(_allocator != 0)
			_allocator->Free(_image);
		_image = 0;
		_allocator = &allocator;
	}
	
	void Load(value_t* image, PolarizationEnum polarization, bool isImaginary)
	{
		std::cout << "Loading " << name(polarization, isImaginary) << '\n';
		if(_imageCount == 1)
			if(_image == 0)
				throw std::runtime_error("Loading image before store");
			else
				memcpy(image, _image, _writer.Width() * _writer.Height() * sizeof(value_t));
		else {
			FitsReader reader(name(polarization, isImaginary));
			reader.Read(image);
		}
	}
	
	void Store(const value_t* image, PolarizationEnum polarization, bool isImaginary)
	{
		std::cout << "Storing " << name(polarization, isImaginary) << '\n';
		if(_imageCount == 1)
		{
			if(_image == 0)
				_image = _allocator->Allocate(_writer.Width() * _writer.Height());
			memcpy(_image, image, _writer.Width() * _writer.Height() * sizeof(value_t));
		}
		else {
			std::string n = name(polarization, isImaginary);
			_writer.Write(n, image);
			_storedNames.insert(n);
		}
	}
	
private:
	std::string name(PolarizationEnum polarization, bool isImaginary) const
	{
		if(isImaginary)
			return _prefix + '-' + Polarization::TypeToShortString(polarization) + "i-tmp.fits";
		else
			return _prefix + '-' + Polarization::TypeToShortString(polarization) + "-tmp.fits";
	}
	FitsWriter _writer;
	size_t _imageCount;
	std::string _prefix;
	
	ImageBufferAllocator<value_t>* _allocator;
	value_t *_image;
	std::set<std::string> _storedNames;
};

#endif
