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
		
		for(std::set<std::string>::const_iterator filenamePtr=_storedNames.begin(); filenamePtr!=_storedNames.end(); ++filenamePtr)
		{
			std::remove(filenamePtr->c_str());
		}
	}
	
	void Initialize(const FitsWriter& writer, size_t polCount, size_t freqCount, const std::string& prefix, ImageBufferAllocator<value_t>& allocator)
	{
		_writer = writer;
		_polCount = polCount;
		_freqCount = freqCount;
		_prefix = prefix;
		if(_allocator != 0)
			_allocator->Free(_image);
		_image = 0;
		_allocator = &allocator;
	}
	
	void SetFitsWriter(const FitsWriter& writer)
	{
		_writer = writer;
	}
	
	void Load(value_t* image, PolarizationEnum polarization, size_t freqIndex, bool isImaginary)
	{
		if(_writer.Width() == 0 || _writer.Height() == 0)
			throw std::runtime_error("Writer is not set.");
		std::cout << "Loading " << name(polarization, freqIndex, isImaginary) << '\n';
		if(_polCount == 1 && _freqCount == 1)
			if(_image == 0)
				throw std::runtime_error("Loading image before store");
			else
				memcpy(image, _image, _writer.Width() * _writer.Height() * sizeof(value_t));
		else {
			FitsReader reader(name(polarization, freqIndex, isImaginary));
			reader.Read(image);
		}
	}
	
	void Store(const value_t* image, PolarizationEnum polarization, size_t freqIndex, bool isImaginary)
	{
		if(_writer.Width() == 0 || _writer.Height() == 0)
			throw std::runtime_error("Writer is not set.");
		std::cout << "Storing " << name(polarization, freqIndex, isImaginary) << '\n';
		if(_polCount == 1 && _freqCount == 1)
		{
			if(_image == 0)
				_image = _allocator->Allocate(_writer.Width() * _writer.Height());
			memcpy(_image, image, _writer.Width() * _writer.Height() * sizeof(value_t));
		}
		else {
			std::string n = name(polarization, freqIndex, isImaginary);
			_writer.Write(n, image);
			_storedNames.insert(n);
		}
	}
	
private:
	std::string name(PolarizationEnum polarization, size_t freqIndex, bool isImaginary) const
	{
		if(_freqCount == 1)
		{
			if(isImaginary)
				return _prefix + '-' + Polarization::TypeToShortString(polarization) + "i-tmp.fits";
			else
				return _prefix + '-' + Polarization::TypeToShortString(polarization) + "-tmp.fits";
		}
		else {
			std::ostringstream str;
			str <<  _prefix + '-' + Polarization::TypeToShortString(polarization);
			if(isImaginary)
				str << 'i';
			str << '-';
			if(freqIndex < 10) str << '0';
			if(freqIndex < 100) str << '0';
			if(freqIndex < 1000) str << '0';
			str << freqIndex << "-tmp.fits";
			return str.str();
		}
	}
	FitsWriter _writer;
	size_t _polCount, _freqCount;
	std::string _prefix;
	
	ImageBufferAllocator<value_t>* _allocator;
	value_t *_image;
	std::set<std::string> _storedNames;
};

#endif
