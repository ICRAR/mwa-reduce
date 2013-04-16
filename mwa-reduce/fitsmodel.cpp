#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "sourcesdf.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: fitsmodel <fitsfile> [spectral index] [ref freq]\n"
			"Prints components in fitsfile to stdout as model.\n";
	} else {
		const char *fitsFilename = argv[1];
		long double spectralIndex, refFreq;
		if(argc > 2)
			spectralIndex = atof(argv[2]);
		else
			spectralIndex = 0.0;
		if(argc > 3)
			refFreq = atof(argv[3]) * 1000000;
		else
			refFreq = 100000000.0;
		FitsReader fitsReader(fitsFilename);
		double *image = new double[fitsReader.ImageWidth() * fitsReader.ImageHeight()];
		fitsReader.Read(image);
		size_t sourceIndex = 0;
		for(size_t y=0; y!=fitsReader.ImageHeight(); ++y)
		{
			for(size_t x=0; x!=fitsReader.ImageWidth(); ++x)
			{
				long double
					l = ((long double) x - (fitsReader.ImageWidth()/2)) * fitsReader.PixelSizeX(),
					m = ((long double) y - (fitsReader.ImageHeight()/2)) * fitsReader.PixelSizeY();
				double value = image[y*fitsReader.ImageWidth() + x];
				if(value != 0.0)
				{
					ModelSource source;
					long double ra, dec;
					ImageCoordinates::LMToRaDec<long double>(l, m, fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), ra, dec);
					std::stringstream nameStr;
					nameStr << "component" << sourceIndex;
					source.SetName(nameStr.str());
					source.SetBrightness(SourceSDFWithSI<long double>(value, spectralIndex, refFreq));
					//std::cout << l << ',' << m << "->" << ra << ',' << dec << '\n';
					source.SetPosRA(ra);
					source.SetPosDec(dec);
					std::cout << source.ToStringLine() << '\n';
					
					++sourceIndex;
				}
			}
		}
		
		delete[] image;
	}
}
