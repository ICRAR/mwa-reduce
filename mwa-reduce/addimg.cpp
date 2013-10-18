#include "fitsreader.h"
#include "fitswriter.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <memory>

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cerr << "Syntax: addimg <outimage> <outweights> <inpimage1> <inpweights1> [<inpimage2> <inpweights2> ...]\n"
			"All images should be fits files. Use -1 as inpweight for unity weight.\n";
	}
	const char *outImageName = argv[1];
	const char *outWeightName = argv[2];
	size_t width = 0, height = 0;
	double *outImage = 0, *outWeights = 0;
	std::unique_ptr<FitsWriter> imgWriter;
	for(int argi=3; argi + 1 < argc; argi += 2)
	{
		const char *inpImageName = argv[argi];
		const char *inpWeightName = argv[argi+1];
		
		FitsReader inpReader(inpImageName);
		if(outImage == 0)
		{
			width = inpReader.ImageWidth(),
			height = inpReader.ImageHeight();
			imgWriter.reset(new FitsWriter(inpReader));
			
			const size_t size = width * height;
			outImage = new double[size];
			outWeights = new double[size];
			for(size_t i=0; i!=size; ++i)
			{
				outImage[i] = 0.0;
				outWeights[i] = 0.0;
			}
		} else {
			if(width != inpReader.ImageWidth() || height != inpReader.ImageHeight())
				throw std::runtime_error("Not all images have same size");
		}
		
		std::vector<double> inpImage(width*height), weightImage(width*height);
		
		inpReader.Read<double>(&inpImage[0]);
		
		if(std::string(inpWeightName) == "-1")
		{
			weightImage.assign(width*height, 1.0);
		}
		else {
			FitsReader weightsReader(inpWeightName);
			if(weightsReader.ImageWidth() != width || weightsReader.ImageHeight() != height)
				throw std::runtime_error("Weights and image do not have same size");
			weightsReader.Read<double>(&weightImage[0]);
		}
			
		// Add the images in
		double *outImagePtr = outImage, *outWeightPtr = outWeights;
		std::vector<double>::iterator inpWeightsIter = weightImage.begin();
		for(std::vector<double>::iterator i=inpImage.begin(); i!=inpImage.end(); ++i)
		{
			*outImagePtr +=  (*i) * (*inpWeightsIter);
			*outWeightPtr += (*inpWeightsIter);
			
			++inpWeightsIter;
			++outImagePtr;
			++outWeightPtr;
		}
	}
	
	// Divide the weight out
	const double *weightsIter = outWeights;
	double *imageEnd = outImage + (width * height);
	for(double *imagePtr=outImage; imagePtr!=imageEnd; ++imagePtr)
	{
		*imagePtr /=  *weightsIter;
		
		++weightsIter;
	}
	
	imgWriter->Write<double>(outImageName, outImage);
	delete[] outImage;
	
	imgWriter->Write<double>(outWeightName, outWeights);
	delete[] outWeights;
}
