#include "fitsreader.h"
#include "fitswriter.h"
#include "matrix2x2.h"
#include "uvector.h"

#include <boost/filesystem/operations.hpp>

#include <stdexcept>
#include <string>
#include <iostream>
#include <complex>
#include <limits>

void read(const FitsReader& templateReader, const std::string& filename, ao::uvector<double>& data)
{
	std::cout << "- " << filename << '\n';
	if(boost::filesystem::exists(filename))
	{
		FitsReader inpReader(filename);
		size_t
			width = inpReader.ImageWidth(),
			height = inpReader.ImageHeight();
			
		if(width != templateReader.ImageWidth() || height != templateReader.ImageHeight())
			throw std::runtime_error("Not all images had same size");
		data.resize(width * height);
		inpReader.Read<double>(&data[0]);
	}
	else {
		std::cout << "Warning: file '" << filename << "' did not exist, assuming it is zero!\n";
		data.assign(templateReader.ImageWidth() * templateReader.ImageHeight(), 0.0);
	}
}

int main(int argc, char *argv[])
{
	if(argc < 5)
	{
		std::cout << "Syntax:\n"
			"  pbaddimg [options] <out-prefix> {<image-prefix> <image-postfix> <beam-prefix> [..]}\n\n"
			"Options:\n"
			"-nopb <nopb-image-prefix> <nopb-weight-prefix>"
			"  Save the images and image weights before beam correction.\n\n"
			"Example:\n"
			"  pbaddimg integrated-stokes wsclean-A image.fits beam-A wsclean-B image.fits beam-B\n"
			"This command will look for image names wsclean-A-XX-image.fits, beam-A-xxr.fits, ... and save to integrated-stokes-i.fits, ...\n";
		return -1;
	}
	
	std::string nopbImagePrefix, nopbWeightPrefix;
	size_t argi = 1;
	while(argv[argi][0] == '-')
	{
		const std::string param(&argv[argi][1]);
		if(param == "nopb")
		{
			++argi;
			nopbImagePrefix = argv[argi];
			++argi;
			nopbWeightPrefix = argv[argi];
		}
		else
			throw std::runtime_error("Invalid command line parameter: "+param);
		++argi;
	}
	
	const std::string outPrefix = argv[argi];
	++argi;
	
	size_t inputCount = (argc - argi) / 3;
	
	ao::uvector<std::complex<double>> outputJonesData, beamJonesData;
	
	std::string firstFilename(std::string(argv[argi])+ "-XX-" + argv[argi+1]);
	FitsReader firstImage(firstFilename);
	const size_t 
		width = firstImage.ImageWidth(),
		height = firstImage.ImageHeight(),
		imgSize = width * height;
		
	for(size_t imageSetIndex=0; imageSetIndex!=inputCount; ++imageSetIndex)
	{
		const std::string imagePrefix = argv[argi], imagePostfix = argv[argi+1], beamPrefix = argv[argi+2];
		argi += 3;
		std::cout << "Processing image set with prefix '" << imagePrefix << "' (" << (imageSetIndex+1) << " / " << inputCount << ")...\n";
	
		std::string
			inpFilenames[4] =
				{ 
					imagePrefix + "-XX-" + imagePostfix, imagePrefix + "-XY-" + imagePostfix,
					imagePrefix + "-XYi-" + imagePostfix, imagePrefix + "-YY-" + imagePostfix
				},
			beamFilename[8] =
				{
					beamPrefix+"-xxr.fits", beamPrefix+"-xxi.fits", beamPrefix+"-xyr.fits", beamPrefix+"-xyi.fits",
					beamPrefix+"-yxr.fits", beamPrefix+"-yxi.fits", beamPrefix+"-yyr.fits", beamPrefix+"-yyi.fits"
				};
		ao::uvector<double> inputData[4], beamData[8];
		for(size_t i=0; i!=4; ++i)
			read(firstImage, inpFilenames[i], inputData[i]);
		for(size_t i=0; i!=8; ++i)
			read(firstImage, beamFilename[i], beamData[i]);
		
		if(outputJonesData.empty())
		{
			outputJonesData.assign(imgSize*4, std::complex<double>(0.0));
			beamJonesData.assign(imgSize*4, std::complex<double>(0.0));
		}
		else if(imgSize*4 != outputJonesData.size())
			throw std::runtime_error("Image sizes did not match!");
		
		ao::uvector<std::complex<double>>::iterator jonesIter = outputJonesData.begin(), beamIter = beamJonesData.begin();
		for(size_t i=0; i!=imgSize; ++i)
		{
			std::complex<double> beamValues[4];
			beamValues[0] = std::complex<double>(beamData[0][i], beamData[1][i]);
			beamValues[1] = std::complex<double>(beamData[2][i], beamData[3][i]);
			beamValues[2] = std::complex<double>(beamData[4][i], beamData[5][i]);
			beamValues[3] = std::complex<double>(beamData[6][i], beamData[7][i]);
			
			std::complex<double> imgValues[4];
			imgValues[0] = inputData[0][i];
			imgValues[1] = std::complex<double>(inputData[1][i], inputData[2][i]);
			imgValues[2] = std::conj(imgValues[1]);
			imgValues[3] = inputData[3][i];
			
			// Calculate Flux += w B* V B  (from: w (B* B) B^-1 V B*^-1 (B* B))
			// 'w' is 1 for now.
			std::complex<double> tempValues[4];
			Matrix2x2::ATimesB(tempValues, beamValues, imgValues);
			Matrix2x2::PlusATimesHermB(&*jonesIter, tempValues, beamValues);
			
			// Calculate Weight += w (B* B) (B* B)
			Matrix2x2::HermATimesB(tempValues, beamValues, beamValues);
			Matrix2x2::HermATimesB(beamValues, tempValues, tempValues);
			Matrix2x2::Add(&*beamIter, beamValues);
			
			jonesIter += 4;
			beamIter += 4;
		}
	}
	
	std::cout << "Dividing by weights and primary beams...\n";
	
	ao::uvector<std::complex<double>>::iterator jonesIter = outputJonesData.begin(), beamIter = beamJonesData.begin();
	ao::uvector<double> stokesOutputImages[4];
	for(size_t p=0; p!=4; ++p)
		stokesOutputImages[p].resize(imgSize);
	for(size_t i=0; i!=imgSize; ++i)
	{
		double stokesMatrix[4];
		
		// Calculate: D sum(W*W)^-1
		if(Matrix2x2::Invert(&*beamIter))
		{
			std::complex<double> correctedLinear[4];
			Matrix2x2::ATimesB(correctedLinear, &*beamIter, &*jonesIter);
			Polarization::LinearToStokes(correctedLinear, stokesMatrix);
		}
		else {
			for(size_t p=0; p!=4; ++p)
				stokesMatrix[p] = std::numeric_limits<double>::quiet_NaN();
		}
		
		for(size_t p=0; p!=4; ++p)
			stokesOutputImages[p][i] = stokesMatrix[p];
		
		jonesIter += 4;
		beamIter += 4;
	}
	
	FitsWriter writer(firstImage);
	
	std::string outFilenames[4] =
		{ outPrefix+"-I.fits", outPrefix+"-Q.fits", outPrefix+"-U.fits", outPrefix+"-V.fits" };
	PolarizationEnum
		stokesPols[4] = { Polarization::StokesI, Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV };

	for(size_t p=0; p!=4; ++p)
	{
		writer.SetPolarization(stokesPols[p]);
		writer.Write<double>(outFilenames[p], &stokesOutputImages[p][0]);
	}
}
