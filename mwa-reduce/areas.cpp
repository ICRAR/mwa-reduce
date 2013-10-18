#include "parser/areaparser.h"

#include "areaset.h"
#include "fitsreader.h"
#include "fitswriter.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cout << "Syntax:\n\tareas <areafile> <reference-fits-file> <output-fits-file>\n";
		return -1;
	}
	
	AreaSet areaSet;
	std::ifstream areaFile(argv[1]);
	AreaParser parser;
	parser.Parse(areaSet, areaFile);
	FitsReader reader(argv[2]);
	const size_t
		width = reader.ImageWidth(),
		height = reader.ImageHeight();
	std::vector<float> image(width * height);
	std::vector<float>::iterator iter=image.begin();
	const size_t SUPERSAMPLE = 1;
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			for(size_t yi=0; yi!=SUPERSAMPLE; ++yi)
			{
				for(size_t xi=0; xi!=SUPERSAMPLE; ++xi)
				{
			
					long double l, m, ra, dec;
					ImageCoordinates::XYToLM<long double>(x*SUPERSAMPLE+xi, y*SUPERSAMPLE+yi, reader.PixelSizeX()/SUPERSAMPLE, reader.PixelSizeY()/SUPERSAMPLE, width*SUPERSAMPLE, height*SUPERSAMPLE, l, m);
					ImageCoordinates::LMToRaDec<long double>(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
					
					std::vector<const SkyArea*> areas;
					areaSet.FindAreas(areas, ra, dec);
					
					*iter += areas.size();
				}
			}
			++iter;
		}
	}
	FitsWriter outputWriter(reader);
	outputWriter.Write(argv[3], &image[0]);
	
	return 0;
}
