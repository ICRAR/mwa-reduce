#include "fitsreader.h"
#include "pngfile.h"
#include "uvector.h"

#include <iostream>

void mapToFire(double value, double& red, double& green, double& blue)
{
	if(value < 0.0)
	{
		red = 0.0;
		green = 0.0;
		blue = 0.0;
	}
	else if(value < 0.25)
	{
		red = value * 2.0;
		green = 0.0;
		blue = 0.0;
	}
	else if(value < 0.5)
	{
		red = value * 2.0;
		green = value * 2.0 - 0.5;
		blue = 0.0;
	}
	else if(value < 0.75)
	{
		red = 1.0;
		green = value * 2.0 - 0.5;
		blue = value * 2.0 - 1.0;
	}
	else if(value < 1.0)
	{
		red = 1.0;
		green = 1.0;
		blue = value * 2.0 - 1.0;
	}
	else
	{
		red = 1.0;
		green = 1.0;
		blue = 1.0;
	}
}

int main(int argc, char* argv[])
{
	if(argc < 6)
	{
		std::cout << "Syntax: fits2png <map> <fitsfile> <min> <max> <pngfile>\n";
		return 0;
	}
	std::string mapName = argv[1], pngFilename(argv[5]);
	FitsReader reader(argv[2]);
	double minVal = atof(argv[3]), maxVal = atof(argv[4]);
	
	const size_t width = reader.ImageWidth(), height = reader.ImageHeight();
	
	std::cout << "Reading...\n";
	
	ao::uvector<double> image(width * height);
	reader.Read(image.data());
	
	
	std::cout << "Writing " << pngFilename << "...\n";
	
	PngFile file(pngFilename, width, height);
	file.BeginWrite();

	ao::uvector<double>::const_iterator imagePtr = image.begin();
	for(unsigned y=0;y<height;++y)
	{
		for(unsigned x=0;x<width;++x)	
		{
			double value = *imagePtr;
			double normValue = (value - minVal) / (maxVal - minVal);
			double r, g, b;
			mapToFire(normValue, r, g, b);
			unsigned rInt = floor(r*256.0), gInt = floor(g*256.0), bInt = floor(b*256.0);
			if(rInt > 255) rInt = 255;
			if(gInt > 255) gInt = 255;
			if(bInt > 255) bInt = 255;
			file.PlotDatapoint(x, height - 1 - y, rInt, gInt, bInt, 255);
			
			++imagePtr;
		}
	}
	file.Close();
}
