#include "banddata.h"
#include "lane.h"
#include "layeredimager.h"
#include "fitswriter.h"
#include "wsinversion.h"

#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[])
{
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0;
	size_t nWLayers = 64;
	bool psf = false;
	
	while(argv[argi][0] == '-')
	{
		const char *param = &argv[argi][1];
		if(strcmp(param, "size") == 0)
		{
			imgWidth = atoi(argv[argi+1]);
			imgHeight = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(strcmp(param, "scale") == 0)
		{
			pixelScale = atof(argv[argi+1]) * M_PI / 180.0;
			++argi;
		}
		else if(strcmp(param, "nwlayers") == 0)
		{
			nWLayers = atoi(argv[argi+1]);
			++argi;
		}
		else if(strcmp(param, "psf") == 0)
		{
			psf = true;
		}
		else {
			throw std::runtime_error("Unknown parameter");
		}
		
		++argi;
	}
	
	const char *msName(argv[argi]);
	const char *fitsfileName(argv[argi+1]);
	
	WSInversion inversionAlgorithm;
	inversionAlgorithm.AddMeasurementSetPath(msName);
	inversionAlgorithm.SetImageWidth(imgWidth);
	inversionAlgorithm.SetImageHeight(imgHeight);
	inversionAlgorithm.SetPixelSizeX(pixelScale);
	inversionAlgorithm.SetPixelSizeY(pixelScale);
	inversionAlgorithm.SetWGridSize(nWLayers);
	inversionAlgorithm.SetDoImagePSF(psf);
	
	inversionAlgorithm.Invert();
	
	std::cout << "Writing image... " << std::flush;
	double
		bandStart = inversionAlgorithm.BandStart(),
		bandEnd = inversionAlgorithm.BandEnd();
		
	FitsWriter writer;
	writer.SetImageDimensions(imgWidth, imgHeight, inversionAlgorithm.PhaseCentreRA(), inversionAlgorithm.PhaseCentreDec(), -pixelScale, pixelScale);
	writer.SetFrequency((bandStart + bandEnd) * 0.5, bandEnd - bandStart);
	writer.SetDate(inversionAlgorithm.StartTime());
	writer.Write(fitsfileName, inversionAlgorithm.ImageResult());
	std::cout << "DONE\n";
}
