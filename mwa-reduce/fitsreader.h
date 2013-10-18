#ifndef FITSREADER_H
#define FITSREADER_H

#include <string>

#include <fitsio.h>

#include "polarizationenum.h"

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
		
		double Frequency() const { return _frequency; }
		double Bandwidth() const { return _bandwidth; }
		
		double DateObs() const { return _dateObs; }
	private:
		float readFloatKey(const char *key);
		double readDoubleKey(const char *key);
		bool readFloatKeyIfExists(const char *key, float &dest);
		bool readDoubleKeyIfExists(const char *key, double &dest);
		std::string readStringKey(const char *key);
		
		std::string _filename;
		fitsfile *_fitsPtr;
		
		void initialize();
		void checkStatus(int status, const std::string &operation=std::string());
		
		size_t _imgWidth, _imgHeight;
		double _phaseCentreRA, _phaseCentreDec;
		double _pixelSizeX, _pixelSizeY;
		double _frequency, _bandwidth, _dateObs;
		PolarizationEnum _polarization;
};

#endif
