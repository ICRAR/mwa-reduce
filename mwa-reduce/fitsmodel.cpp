#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "sourcesdf.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "cleanalgorithm.h"

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
		
		Model model;
		CleanAlgorithm::GetModelFromImage(model, image, fitsReader.ImageWidth(), fitsReader.ImageHeight(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), fitsReader.PixelSizeX(), fitsReader.PixelSizeY(), spectralIndex, refFreq);
		
		for(Model::const_iterator i=model.begin(); i!=model.end(); ++i)
		{
			std::cout << i->ToString() << '\n';
		}
		
		delete[] image;
	}
}
