#ifndef FITSWRITER_H
#define FITSWRITER_H

#include <string>

class FitsWriter
{
	public:
		FitsWriter(const std::string &filename) 
		: _filename(filename) { }
		template<typename NumType> void Write(const NumType *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs);
	private:
		std::string _filename;
		
		void checkStatus(int status);
		void julianDateToYMD(double jd, int &year, int &month, int &day);
};

#endif
