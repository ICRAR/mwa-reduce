#ifndef FITSWRITER_H
#define FITSWRITER_H

#include <string>

class FitsWriter
{
	public:
		FitsWriter(const std::string &filename) 
		: _filename(filename) { }
		template<typename NumType> void Write(NumType *image, size_t imgSize, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY);
	private:
		std::string _filename;
		
		void checkStatus(int status);
};

#endif
