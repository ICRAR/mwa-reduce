#include "cleanalgorithm.h"
#include "inversionalgorithm.h"
#include "wsinversion.h"
#include "fitswriter.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << "Syntax:\tclean [options] <input-ms> <image-prefix>\n"
			"Will create cleaned images of the input ms. DATA column will be used by default.\n"
			"Options can be:\n"
			"\t-size <width> <height>\n"
			"\t-scale <pixel-scale>\n"
			"\t   Scale of a pixel in degrees, e.g. 0.012.\n"
			"\t-nwlayers <nwlayers>\n"
			"\t   Number of w-layers to use\n"
			"\t-niter <niter>\n"
			"\t   Maximum number of clean iterations to perform\n"
			"\t-threshold <threshold>\n"
			"\t   Stopping clean thresholding in Jy\n"
			"\t-gain <gain>\n"
			"\t   Cleaning gain: Ratio of peak that will be subtracted in each iteration\n"
			"\t-pol <xx, yy, xy, yx or stokesi>\n";
		return -1;
	}
	
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0, threshold = 0.0, gain = 0.1;
	size_t nWLayers = 64, nIter = 500;
	enum InversionAlgorithm::PolarizationEnum polarization = InversionAlgorithm::StokesI;
	
	while(argv[argi][0] == '-')
	{
		const char *param = &argv[argi][1];
		if(strcmp(param, "size") == 0)
		{
			imgWidth = atoi(argv[argi+1]);
			imgHeight = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(strcmp(param, "scale") == 0)
		{
			pixelScale = atof(argv[argi+1]) * M_PI / 180.0;
			++argi;
		}
		else if(strcmp(param, "nwlayers") == 0)
		{
			nWLayers = atoi(argv[argi+1]);
			++argi;
		}
		else if(strcmp(param, "gain") == 0)
		{
			gain = atof(argv[argi + 1]);
			++argi;
		}
		else if(strcmp(param, "niter") == 0)
		{
			nIter = atoi(argv[argi + 1]);
			++argi;
		}
		else if(strcmp(param, "threshold") == 0)
		{
			threshold = atof(argv[argi + 1]);
			++argi;
		}
		else if(strcmp(param, "pol") == 0)
		{
			std::string polStr = argv[argi + 1];
			boost::to_lower(polStr);
			if(polStr == "xx")
				polarization = InversionAlgorithm::XX;
			else if(polStr == "xy")
				polarization = InversionAlgorithm::XY;
			else if(polStr == "yx")
				polarization = InversionAlgorithm::YX;
			else if(polStr == "yy")
				polarization = InversionAlgorithm::YY;
			else if(polStr == "stokesi")
				polarization = InversionAlgorithm::StokesI;
			++argi;
		}
		else {
			throw std::runtime_error("Unknown parameter");
		}
		
		++argi;
	}
	
	const char *msName(argv[argi]);
	const char *fileNamePrefix(argv[argi+1]);
	
	std::unique_ptr<InversionAlgorithm> inversionAlgorithm(new WSInversion());
	inversionAlgorithm->SetMeasurementSetPath(msName);
	inversionAlgorithm->SetImageWidth(imgWidth);
	inversionAlgorithm->SetImageHeight(imgHeight);
	inversionAlgorithm->SetPixelSizeX(pixelScale);
	inversionAlgorithm->SetPixelSizeY(pixelScale);
	inversionAlgorithm->SetWGridSize(nWLayers);
	inversionAlgorithm->SetPolarization(polarization);
	
	inversionAlgorithm->SetDoImagePSF(true);
	inversionAlgorithm->Execute();
	std::vector<double> psf(imgWidth * imgHeight);
	memcpy(&psf[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
	const double
		ra = inversionAlgorithm->ImageResultRA(),
		dec = inversionAlgorithm->ImageResultDec();
	
	std::cout << "Writing psf image... " << std::flush;
	FitsWriter psfWriter(std::string(fileNamePrefix) + "-psf.fits");
	psfWriter.Write(&psf[0], imgWidth, imgHeight, ra, dec, -pixelScale, pixelScale);
	std::cout << "DONE\n";
	
	inversionAlgorithm->SetDoImagePSF(false);
	inversionAlgorithm->Execute();
	std::vector<double> model(imgWidth * imgHeight), residual(imgWidth * imgHeight);
	memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
	inversionAlgorithm.reset();
	
	CleanAlgorithm cleanAlgorithm;
	cleanAlgorithm.SetMaxNIter(nIter);
	cleanAlgorithm.SetThreshold(threshold);
	cleanAlgorithm.SetSubtractionGain(gain);
	cleanAlgorithm.ExecuteMajorIteration(&residual[0], &model[0], &psf[0], imgWidth, imgHeight);
	
	std::cout << "Writing residual image... " << std::flush;
	FitsWriter resWriter(std::string(fileNamePrefix) + "-residual.fits");
	resWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, -pixelScale, pixelScale);
	std::cout << "DONE\n";
	
	std::cout << "Writing model image... " << std::flush;
	FitsWriter modelWriter(std::string(fileNamePrefix) + "-model.fits");
	modelWriter.Write(&model[0], imgWidth, imgHeight, ra, dec, -pixelScale, pixelScale);
	std::cout << "DONE\n";
}
