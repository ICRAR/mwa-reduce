#include "fitsreader.h"
#include "fitswriter.h"
#include "uvector.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <memory>

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cerr << "Syntax: addimg <outimage> <outweights> <inpimage1> <inpbeam1> [<inpimage2> <inpbeam2> ...]\n"
			"All images should be fits files. Use -1 as inpweight for unity weight, or \"-c inpbeam-real inpbeam-imag\" for complex beam.\n";
	}
	const char *outImageName = argv[1];
	const char *outWeightName = argv[2];
	size_t width = 0, height = 0, count = 0;
	double *outImage = 0, *outWeights = 0;
	std::unique_ptr<FitsWriter> imgWriter;
	double frequencySum = 0.0, lowestFreq = 0.0, highestFreq = 0.0;
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
		
		double thisFrequency = inpReader.Frequency();
		frequencySum += thisFrequency;
		if(count == 0)
		{
			lowestFreq = thisFrequency - inpReader.Bandwidth()*0.5;
			highestFreq = thisFrequency + inpReader.Bandwidth()*0.5;
		}
		else {
			if(thisFrequency - inpReader.Bandwidth()*0.5 < lowestFreq)
				lowestFreq = thisFrequency - inpReader.Bandwidth()*0.5;
			if(thisFrequency + inpReader.Bandwidth()*0.5 > highestFreq)
				highestFreq = thisFrequency + inpReader.Bandwidth()*0.5;
		}
		count++;
		
		ao::uvector<double> inpImage(width*height), weightImage(width*height);
		
		inpReader.Read<double>(&inpImage[0]);
		
		if(std::string(inpWeightName) == "-1")
		{
			weightImage.assign(width*height, 1.0);
		}
		else if(std::string(inpWeightName) == "-c")
		{
			argi += 2;
			FitsReader realReader(argv[argi]), imagReader(argv[argi+1]);
			if(realReader.ImageWidth() != width || realReader.ImageHeight() != height)
				throw std::runtime_error("Real beam and image do not have same size");
			if(imagReader.ImageWidth() != width || imagReader.ImageHeight() != height)
				throw std::runtime_error("Imaginary beam and image do not have same size");
			
			ao::uvector<double> realImage(width*height), imagImage(width*height);
			realReader.Read<double>(&realImage[0]);
			imagReader.Read<double>(&imagImage[0]);
			for(size_t j=0; j!=width*height; ++j)
			{
				double r = realImage[j], i = imagImage[j];
				weightImage[j] = r*r + i*i;
			}
		}
		else {
			FitsReader weightsReader(inpWeightName);
			if(weightsReader.ImageWidth() != width || weightsReader.ImageHeight() != height)
				throw std::runtime_error("Weights and image do not have same size");
			weightsReader.Read<double>(&weightImage[0]);
		}
			
		// Add the images in
		double *outImagePtr = outImage, *outWeightPtr = outWeights;
		ao::uvector<double>::iterator inpWeightsIter = weightImage.begin();
		for(ao::uvector<double>::iterator i=inpImage.begin(); i!=inpImage.end(); ++i)
		{
			double beamVal = *inpWeightsIter;
			*outImagePtr +=  (*i) * beamVal;
			*outWeightPtr += beamVal * beamVal;
			
			++inpWeightsIter;
			++outImagePtr;
			++outWeightPtr;
		}
		
		std::cout << '.' << std::flush;
	}
	std::cout << '\n';
	
	// Divide the weight out
	const double *weightsIter = outWeights;
	double *imageEnd = outImage + (width * height);
	for(double *imagePtr=outImage; imagePtr!=imageEnd; ++imagePtr)
	{
		double weight = *weightsIter;
		*imagePtr /=  weight;
		
		++weightsIter;
	}
	
	imgWriter->SetFrequency(frequencySum / count, (highestFreq - lowestFreq));
	imgWriter->Write<double>(outImageName, outImage);
	delete[] outImage;
	
	imgWriter->Write<double>(outWeightName, outWeights);
	delete[] outWeights;
}
