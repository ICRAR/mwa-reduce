#ifndef FITSREADER_H
#define FITSREADER_H

#include <string>

#include <fitsio.h>

class FitsReader
{
	public:
		FitsReader(const std::string &filename) 
		: _filename(filename) { initialize(); }
		~FitsReader();
		
		template<typename NumType> void Read(NumType *image);
		
		size_t ImageWidth() const { return _imgWidth; }
		size_t ImageHeight() const { return _imgHeight; }
		
		double PhaseCentreRA() const { return _phaseCentreRA; }
		double PhaseCentreDec() const { return _phaseCentreDec; }
		
		double PixelSizeX() const { return _pixelSizeX; }
		double PixelSizeY() const { return _pixelSizeY; }
	private:
		float readFloatKey(const char *key);
		void readFloatKeyIfExists(const char *key, float &dest);
		std::string readStringKey(const char *key);
		
		std::string _filename;
		fitsfile *_fitsPtr;
		
		void initialize();
		void checkStatus(int status, const std::string &operation=std::string());
		
		size_t _imgWidth, _imgHeight;
		double _phaseCentreRA, _phaseCentreDec;
		double _pixelSizeX, _pixelSizeY;
};

#endif
