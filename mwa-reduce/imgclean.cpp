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
		std::cerr << "Syntax: imgclean <inpimage> [<psfimage> <outpimage> <modelimage>]\n"
			"All images should be fits files. If only an inpimage is specified, the first peak is returned.\n";
		return -1;
	}
	bool onlyFindPeak = argc == 2;
	const char *inpImageName = argv[1];
	const char *psfImageName, *outImageName, *modelImageName;
	if(onlyFindPeak)
	{
		psfImageName = 0;
		outImageName = 0;
		modelImageName = 0;
	} else {
		psfImageName = argv[2];
		outImageName = argv[3];
		modelImageName = argv[4];
	}
	double threshold = 5.0, subtractionFactor = 0.25;
	size_t maxIter = 50000;
	bool allowNegativeComponents = true, stopWhenDiverging = false;
	
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
	double peak = SimpleClean::FindPeak(&image[0], width, height, componentX, componentY, allowNegativeComponents, 0, height, 0.05);
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
		double lastPeak = peak;
		std::cout << "Initial peak: " << peak << '\n';
		size_t iterationNumber = 0;
		FitsReader psfReader(psfImageName);
		psfReader.Read<double>(&psf[0]);
		while(fabs(peak) > threshold && iterationNumber < maxIter && (!stopWhenDiverging || peak<=lastPeak))
		{
			if(iterationNumber % 10 == 0)
				std::cout << "Iteration " << iterationNumber << ": (" << componentX << ',' << componentY << "), " << peak << " Jy\n";
			SimpleClean::SubtractImage(&image[0], &psf[0], width, height, componentX, componentY, subtractionFactor * peak);
			model[componentX + componentY*width] += subtractionFactor * peak;
			
			lastPeak = peak;
			peak = SimpleClean::FindPeak(&image[0], width, height, componentX, componentY, allowNegativeComponents, 0, height, 0.05);
			++iterationNumber;
		}
		std::cout << "Stopped on peak " << peak << '\n';
		
		FitsWriter imgWriter(inpReader);
		imgWriter.Write<double>(outImageName, &image[0]);
		imgWriter.Write<double>(modelImageName, &model[0]);
	}
}
