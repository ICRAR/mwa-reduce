#include <iostream>
#include <stdexcept>
#include <vector>
#include "fitsreader.h"
#include "model.h"
#include "modelrenderer.h"
#include "fitswriter.h"

int main(int argc, char* argv[])
{
	if(argc == 1)
		std::cout << "syntax: render [-t templatefits] [-b] <model> <outputfits>\n";
	else {
		std::string templateFits;
		
		int argi = 1;
		while(argi < argc && argv[argi][0] == '-')
		{
			std::string param(&argv[argi][1]);
			if(param == "t") {
				++argi;
				templateFits = argv[argi];
			}
			else throw std::runtime_error("Invalid param");
			++argi;
		}
	
		Model model(argv[argi]);
		std::string outputFitsName(argv[argi+1]);
	
		size_t width = 1024, height = 1024;
		double ra = 0.0, dec = 0.0;
		double pixelSizeX = 0.01, pixelSizeY = 0.01;
		double bandwidth = 1000000.0, dateObs = 0.0, frequency = 150000000.0;
		
		if(!templateFits.empty())
		{
			FitsReader reader(templateFits);
			width = reader.ImageWidth();
			height = reader.ImageHeight();
			ra = reader.PhaseCentreRA();
			dec = reader.PhaseCentreDec();
			pixelSizeX = reader.PixelSizeX();
			pixelSizeY = reader.PixelSizeY();
			bandwidth = reader.Bandwidth();
			dateObs = reader.DateObs();
			frequency = reader.Frequency();
		}
		
		std::vector<double> image(width * height);
		ModelRenderer renderer(ra, dec, pixelSizeX, pixelSizeY);
		renderer.RenderModel(&image[0], width, height, model, frequency-bandwidth*0.5, frequency+bandwidth*0.5, 0);
		
		FitsWriter writer(outputFitsName);
		writer.Write(&image[0], width, height, ra, dec, pixelSizeX, pixelSizeY, frequency, bandwidth, dateObs);
	}
}
