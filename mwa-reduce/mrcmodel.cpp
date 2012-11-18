#include <iostream>
#include <cmath>
#include <fstream>

#include "sourcestrength.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "mrccatalogue.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: mrcmodel <mrc-catalogue> <fitsfile> [<flux threshold> <ref freq>]\n"
			"Prints model to stdout, constructed from given mrc catalogue and\n"
			"image fitsfile.\n";
	} else {
		const char *catFilename = argv[1];
		const char *fitsFilename = argv[2];
		long double threshold;
		if(argc > 3)
			threshold = atof(argv[3]);
		else
			threshold = 0.0;
		long double refFreq;
		if(argc > 4)
			refFreq = atof(argv[4]);
		else
			refFreq = 408.0;
		
		MRCCatalogue catalogue(catFilename);
		FitsReader fitsReader(fitsFilename);
		double *image = new double[fitsReader.ImageWidth() * fitsReader.ImageHeight()];
		fitsReader.Read(image);
		
		ModelSource source;
		while(catalogue.ReadNext(source))
		{
			long double l, m;
			ImageCoordinates::RaDecToLM<long double>(source.PosRA(), source.PosDec(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), l, m);
			double
				x = l / fitsReader.PixelSizeX() + fitsReader.ImageWidth()/2,
				y = m / fitsReader.PixelSizeY() + fitsReader.ImageHeight()/2;
			//std::cout << x << ',' << y << '\n';
			bool aboveThreshold = (source.Brightness().FluxAtFrequency(refFreq) > threshold);
			if(aboveThreshold && x > 0.0 && y > 0.0 && x < fitsReader.ImageWidth() && y < fitsReader.ImageHeight())
			{
				double value = 0.0;
				bool sourceStrengthFromFits = false;
				if(sourceStrengthFromFits)
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
					
					source.Brightness() = SourceStrength<long double>(value, 0.0, 1.0);
				}
				
				std::cout << source.ToStringLine() << '\n';
			}
		}
		
		delete[] image;
		
	}
}
