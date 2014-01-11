#include <iostream>
#include "fitsreader.h"
#include "model.h"
#include "ioninterpolator.h"
#include "ionsolutionfile.h"

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		std::cout << "Syntax:\n\tapplyion <input fits> <output fits> <model> <ion-solutions>\n";
	}
	const char
		*inputFilename = argv[1],
		*outputFilename = argv[2],
		*modelFilename = argv[3],
		*solutionsFilename = argv[4];
		
	FitsReader reader(inputFilename);
	Model model(modelFilename);
	IonInterpolator interpolator(model, reader);
	IonSolutionFile solutions(solutionsFilename);
	interpolator.Initialize(solutions, 0, solutions.IntervalCount(), 0, solutions.ChannelCount(), 0);
	
	const size_t width = reader.ImageWidth(), height = reader.ImageHeight();
	ao::uvector<double>
		gainImage(width * height),
		dlImage(width * height),
		dmImage(width * height),
		outImage(width * height);
	interpolator.Interpolate(gainImage.data(), solutions, IonSolutionFile::GainSolution);
	interpolator.Interpolate(dlImage.data(), solutions, IonSolutionFile::DlSolution);
	interpolator.Interpolate(dmImage.data(), solutions, IonSolutionFile::DmSolution);
	
	double
		*gainPtr = gainImage.data(),
		*dlPtr = dlImage.data(),
		*dmPtr = dmImage.data(),
		*outPtr = outImage.data();
	for(size_t y=0; y!=height; ++y)
	{
		for(size_t x=0; x!=width; ++x)
		{
			const double
				dl = *dlPtr,
				dm = *dmPtr,
				gain = *gainPtr;
			
			double l, m;
			ImageCoordinates::XYToLM(x, y, reader.PixelSizeX(), reader.PixelSizeY(), width, height, l, m);
			
			
			++gainPtr;
			++dlPtr;
			++dmPtr;
			++outPtr;
		}
	}
}
