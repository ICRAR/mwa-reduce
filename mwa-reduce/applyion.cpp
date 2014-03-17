#include <iostream>
#include <limits>

#include "fitsreader.h"
#include "fitswriter.h"
#include "ioninterpolator.h"
#include "ionsolutionfile.h"
#include "model.h"

double sample(const double* image, size_t width, size_t height, double x, double y)
{
	double xleft = floor(x), ytop = floor(y);
	if(xleft < 0.0 || ytop < 0.0) return std::numeric_limits<double>::quiet_NaN();
	
	size_t x1 = size_t(xleft), y1 = size_t(ytop);
	if(x1+1>=width || y1+1>=height) return std::numeric_limits<double>::quiet_NaN();
	
	size_t index = x1 + y1*width;
	
	// Do bilinair interpolation
	double x1dist = x - xleft, y1dist = y - ytop, x2dist = 1.0 - x1dist;
	double
		v11 = image[index],
		v12 = image[index+1],
		v21 = image[index+width],
		v22 = image[index+width+1];
	return
		(v11 * x2dist + v12 * x1dist) * (1.0 - y1dist) +
		(v21 * x2dist + v22 * x1dist) * y1dist;
}

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		std::cout << "Syntax:\n\tapplyion <input fits> <output fits> <model> <ion-solutions>\n";
		return -1;
	}
	const char
		*inputFilename = argv[1],
		*outputFilename = argv[2],
		*modelFilename = argv[3],
		*solutionsFilename = argv[4];
		
	FitsReader reader(inputFilename);
	Model model(modelFilename);
	IonInterpolator interpolator(model, reader);
	IonSolutionFile solutions;
	solutions.OpenForReading(solutionsFilename);
	std::cout << "Model has " << model.ComponentCount() << " components.\n";
	std::cout << "Solutions have " << solutions.DirectionCount() << " directions.\n";
	if(model.SourceCount() != solutions.DirectionCount())
		throw std::runtime_error("Nr of clusters in model does not match number of solution directions!");
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
	
	ao::uvector<double> image(width * height);
	reader.Read(image.data());
	
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
			l += dl;
			m += dm;
			
			double xf, yf;
			ImageCoordinates::LMToXYfloat(l, m, reader.PixelSizeX(), reader.PixelSizeY(), width, height, xf, yf);
			
			double val = sample(image.data(), width, height, xf, yf);
			*outPtr = val / gain;
			
			++gainPtr;
			++dlPtr;
			++dmPtr;
			++outPtr;
		}
	}
	
	FitsWriter writer(reader);
	writer.Write(outputFilename, outImage.data());
}
