#include <iostream>
#include <cmath>
#include <fstream>

#include "sourcestrength.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"

int main(int argc, char **argv)
{
	if(argc != 3)
	{
		std::cout << "Usage: apparently <model> <fitsfile>\n"
			"Prints apparent model to stdout, constructed from given model and\n"
			"image fitsfile.\n";
	} else {
		const char *modelFilename = argv[1];
		const char *fitsFilename = argv[2];
		
		Model model(modelFilename);
		FitsReader fitsReader(fitsFilename);
		double *image = new double[fitsReader.ImageWidth() * fitsReader.ImageHeight()];
		fitsReader.Read(image);
		
		for(Model::const_iterator s=model.begin();s!=model.end();++s)
		{
			ModelSource source = *s;
			long double l, m;
			ImageCoordinates::RaDecToLM<long double>(source.PosRA(), source.PosDec(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), l, m);
			double
				x = l / fitsReader.PixelSizeX() + fitsReader.ImageWidth()/2,
				y = m / fitsReader.PixelSizeY() + fitsReader.ImageHeight()/2;
			double value = 0.0;
			//std::cout << x << ',' << y << '\n';
			if(x > 0.0 && y > 0.0 && x < fitsReader.ImageWidth() && y < fitsReader.ImageHeight())
			{
				size_t xi = size_t(floor(x)), yi = size_t(floor(y));
				value = image[yi*fitsReader.ImageWidth() + xi];
				
				++xi;
				if(xi < fitsReader.ImageWidth()) {
					value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
					
					++yi;
					if(yi < fitsReader.ImageHeight())
						value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
				}
				
				--xi;
				if(yi < fitsReader.ImageHeight())
					value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
			}
			if(value > 0.0) {
				source.Brightness() = SourceStrength<long double>(value, 0.0, 1.0);
				std::cout << source.ToStringLine() << '\n';
			} else {
				std::cout << '#' << source.ToStringLine() << '\n';
			}
		}
		
		delete[] image;
		
	}
}
