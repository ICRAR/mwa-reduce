#include "cleanalgorithms/simpleclean.h"

#include "fitsreader.h"
#include "fitswriter.h"
#include "imagecoordinates.h"
#include "modelsource.h"
#include "model.h"

#include <vector>
#include <stdexcept>

#include <iostream>

int main(int argc, char *argv[])
{
	if(argc != 5 && argc != 2)
	{
		std::cerr << "Syntax: imgclean [options] <inpimage> [<psfimage> <outpimage> <modelimage>]\n"
			"All images should be fits files. If only an inpimage is specified, the first peak is returned.\n"
			"Options:\n"
			"-t <threshold>\n"
			"  Stop at given threshold\n"
			"-g <gain>\n"
			"  Use given loop gain\n"
			"-b <border size>\n"
			"  Specify border size that won't be cleaned, as ratio of image size (after boxing).\n"
			"-box <x> <y> <width> <height>\n"
			"  Cut the inpimage centred on x,y with given width and height and only clean that.\n";
		return -1;
	}
	int argi=1;
	double threshold = 1.0, gain = 0.1, border = 0.05;
	size_t boxX = 0, boxY = 0, boxWidth = 0, boxHeight = 0;
	while(argv[argi][0] == '-')
	{
		const std::string param = &argv[argi][1];
		if(param == "box")
		{
		}
		else if(param == "b") {
			++argi;
			border = atof(argv[argi]);
		}
		else if(param == "t") {
			++argi;
			threshold = atof(argv[argi]);
		}
		else if(param == "g") {
			++argi;
			gain = atof(argv[argi]);
		}
		else
			throw std::runtime_error("Unknown parameter");
		++argi;
	}
	
	bool onlyFindPeak = (argc-argi) == 1;
	const char *inpImageName = argv[argi];
	const char *psfImageName, *outImageName, *modelImageName;
	if(onlyFindPeak)
	{
		psfImageName = 0;
		outImageName = 0;
		modelImageName = 0;
	} else {
		psfImageName = argv[argi+1];
		outImageName = argv[argi+2];
		modelImageName = argv[argi+3];
	}
	size_t maxIter = 50000;
	bool allowNegativeComponents = true;
	
	FitsReader inpReader(inpImageName);
	const size_t width = inpReader.ImageWidth();
	const size_t height = inpReader.ImageHeight();
		
	const size_t size = width * height;
	std::vector<double>
		image(size),
		psf(size),
		model(size, 0.0);
	
	inpReader.Read<double>(&image[0]);
	
	size_t componentX, componentY;
	double peak = SimpleClean::FindPeak(&image[0], width, height, componentX, componentY, allowNegativeComponents, 0, height, border);
	if(onlyFindPeak)
	{
		double l, m, ra, dec;
		ImageCoordinates::XYToLM(componentX, componentY, inpReader.PixelSizeX(), inpReader.PixelSizeY(), width, height, l, m);
		ImageCoordinates::LMToRaDec(l, m, inpReader.PhaseCentreRA(), inpReader.PhaseCentreDec(), ra, dec);
		ModelComponent component;
		component.SetPosRA(ra);
		component.SetPosDec(dec);
		component.SetSED(SpectralEnergyDistribution(peak, 1.0));
		ModelSource source;
		source.SetName("clcomp");
		source.AddComponent(component);
		Model outputModel;
		outputModel.AddSource(source);
		outputModel.Save(std::cout);
	} else {
		std::cout << "Initial peak: " << peak << '\n';
		size_t iterationNumber = 0;
		FitsReader psfReader(psfImageName);
		psfReader.Read<double>(&psf[0]);
		while(fabs(peak) > threshold && iterationNumber < maxIter)
		{
			if(iterationNumber % 10 == 0)
				std::cout << "Iteration " << iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
			SimpleClean::SubtractImage(&image[0], &psf[0], width, height, componentX, componentY, gain * peak);
			model[componentX + componentY*width] += gain * peak;
			
			peak = SimpleClean::FindPeak(&image[0], width, height, componentX, componentY, allowNegativeComponents, 0, height, border);
			++iterationNumber;
		}
		std::cout << "Stopped on peak " << peak << '\n';
		
		FitsWriter imgWriter(inpReader);
		imgWriter.Write<double>(outImageName, &image[0]);
		imgWriter.Write<double>(modelImageName, &model[0]);
	}
}
