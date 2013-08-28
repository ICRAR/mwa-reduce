#include "cleanalgorithm.h"
#include "inversionalgorithm.h"
#include "wsinversion.h"
#include "fitswriter.h"
#include "modelrenderer.h"
#include "model.h"
#include "areaset.h"
#include "parser/areaparser.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout << "Syntax:\twsclean [options] <input-ms> [<2nd-ms> [..]]\n"
			"Will create cleaned images of the input ms(es). DATA column will be used by default.\n"
			"If multiple mses are specified, they need to be phase-rotated to the same point on the sky.\n"
			"Options can be:\n"
			"\t-name <image-prefix>\n"
			"\t   Use image-prefix as prefix for output files. Default is 'wsclean'.\n"
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
			"\t   Cleaning gain: Ratio of peak that will be subtracted in each iteration (default = 0.1).\n"
			"\t-mgain <gain>\n"
			"\t   Cleaning gain for major iterations: Ratio of peak that will be subtracted in each major\n"
			"\t   iteration (default = 1.0, to use major iterations, 0.9 is a good value).\n"
			"\t-pol <xx, yy, xy, yx or stokesi>\n"
			"\t-negative\n"
			"\t   Allow negative components during cleaning\n"
			"\t-datacolumn <columnname>\n"
			"\t-addmodel <modelfile>\n"
			"\t-savemodel <modelfile>\n";
		return -1;
	}
	
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0, threshold = 0.0, gain = 0.1, mGain = 1.0;
	size_t nWLayers = 64, nIter = 500;
	std::string columnName = "DATA", addModelFilename, saveModelFilename, cleanAreasFilename;
	enum InversionAlgorithm::PolarizationEnum polarization = InversionAlgorithm::StokesI;
	std::string prefixName = "wsclean";
	bool allowNegative = false;
	enum LayeredImager::GridModeEnum gridMode = LayeredImager::NearestNeighbour;
	
	while(argi < argc && argv[argi][0] == '-')
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
			++argi;
			pixelScale = atof(argv[argi]) * M_PI / 180.0;
		}
		else if(strcmp(param, "nwlayers") == 0)
		{
			++argi;
			nWLayers = atoi(argv[argi]);
		}
		else if(strcmp(param, "gain") == 0)
		{
			++argi;
			gain = atof(argv[argi]);
		}
		else if(strcmp(param, "mgain") == 0)
		{
			++argi;
			mGain = atof(argv[argi]);
		}
		else if(strcmp(param, "niter") == 0)
		{
			++argi;
			nIter = atoi(argv[argi]);
		}
		else if(strcmp(param, "threshold") == 0)
		{
			++argi;
			threshold = atof(argv[argi]);
		}
		else if(strcmp(param, "datacolumn") == 0)
		{
			++argi;
			columnName = argv[argi];
		}
		else if(strcmp(param, "pol") == 0)
		{
			++argi;
			std::string polStr = argv[argi];
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
		}
		else if(strcmp(param, "negative") == 0)
		{
			allowNegative = true;
		}
		else if(strcmp(param, "addmodel") == 0)
		{
			++argi;
			addModelFilename = argv[argi];
		}
		else if(strcmp(param, "savemodel") == 0)
		{
			++argi;
			saveModelFilename = argv[argi];
		}
		else if(strcmp(param, "cleanareas") == 0)
		{
			++argi;
			cleanAreasFilename = argv[argi];
		}
		else if(strcmp(param, "name") == 0)
		{
			++argi;
			prefixName = argv[argi];
		}
		else if(strcmp(param, "gridmode") == 0)
		{
			++argi;
			std::string gridModeStr = argv[argi];
			boost::to_lower(gridModeStr);
			if(gridModeStr == "kb" || gridModeStr == "kaiserbessel" || gridModeStr == "kaiser-bessel")
				gridMode = LayeredImager::KaiserBessel;
			else if(gridModeStr == "nn" || gridModeStr == "nearestneighbour")
				gridMode = LayeredImager::NearestNeighbour;
			else
				throw std::runtime_error("Invalid gridding mode: should be either kb (Kaiser-Bessel) or nn (NearestNeighbour)");
		}
		else {
			throw std::runtime_error("Unknown parameter");
		}
		
		++argi;
	}
	
	if(argi == argc)
		throw std::runtime_error("No input measurement sets given.");
	
	std::unique_ptr<InversionAlgorithm> inversionAlgorithm(new WSInversion());
	static_cast<WSInversion&>(*inversionAlgorithm).SetGridMode(gridMode);
	
	for(int i=argi; i != argc; ++i)
		inversionAlgorithm->AddMeasurementSetPath(argv[i]);
	
	inversionAlgorithm->SetImageWidth(imgWidth);
	inversionAlgorithm->SetImageHeight(imgHeight);
	inversionAlgorithm->SetPixelSizeX(pixelScale);
	inversionAlgorithm->SetPixelSizeY(pixelScale);
	inversionAlgorithm->SetWGridSize(nWLayers);
	inversionAlgorithm->SetPolarization(polarization);
	inversionAlgorithm->SetDataColumnName(columnName);
	
	double
		ra, dec,
		freqHigh, freqLow, freqCentre,
		bandwidth, beamSize, dateObs;
		
	std::vector<double> psf;
	
	if(nIter > 0)
	{
		std::cout << " == Constructing PSF ==\n";
		inversionAlgorithm->SetDoImagePSF(true);
		inversionAlgorithm->Invert();
		psf.resize(imgWidth * imgHeight);
		memcpy(&psf[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
		
		ra = inversionAlgorithm->ImageResultRA(),
		dec = inversionAlgorithm->ImageResultDec(),
		freqHigh = inversionAlgorithm->ImageHighestFrequencyChannel(),
		freqLow = inversionAlgorithm->ImageLowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = inversionAlgorithm->ImageBandEnd() - inversionAlgorithm->ImageBandStart(),
		beamSize = inversionAlgorithm->ImageBeamSize(),
		dateObs = inversionAlgorithm->ImageStartTime();
		
		std::cout << "Writing psf image... " << std::flush;
		FitsWriter psfWriter(std::string(prefixName) + "-psf.fits");
		psfWriter.Write(&psf[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
		std::cout << "DONE\n";
	
		if(inversionAlgorithm->HasGriddingCorrectionImage())
		{
			std::cout << "Writing gridding correction image... " << std::flush;
			FitsWriter griddingWriter(std::string(prefixName) + "-gridding.fits");
			std::vector<double> gridding(imgWidth * imgHeight);
			inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
			griddingWriter.Write(&gridding[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
			std::cout << "DONE\n";
		}
	}
	
	std::cout << " == Constructing image ==\n";
	inversionAlgorithm->SetDoImagePSF(false);
	inversionAlgorithm->Invert();
	
	ra = inversionAlgorithm->ImageResultRA(),
	dec = inversionAlgorithm->ImageResultDec(),
	freqHigh = inversionAlgorithm->ImageHighestFrequencyChannel(),
	freqLow = inversionAlgorithm->ImageLowestFrequencyChannel(),
	freqCentre = (freqHigh + freqLow) * 0.5,
	bandwidth = inversionAlgorithm->ImageBandEnd() - inversionAlgorithm->ImageBandStart(),
	beamSize = inversionAlgorithm->ImageBeamSize(),
	dateObs = inversionAlgorithm->ImageStartTime();
	
	std::vector<double> modelImage(imgWidth * imgHeight), residual(imgWidth * imgHeight);
	memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
	
	if(mGain == 1.0)
		inversionAlgorithm.reset();
	
	std::cout << "Writing dirty image... " << std::flush;
	FitsWriter dirtyWriter(std::string(prefixName) + "-dirty.fits");
	dirtyWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
	std::cout << "DONE\n";
	
	if(nIter > 0)
	{
		CleanAlgorithm cleanAlgorithm;
		cleanAlgorithm.SetMaxNIter(nIter);
		cleanAlgorithm.SetThreshold(threshold);
		cleanAlgorithm.SetSubtractionGain(gain);
		cleanAlgorithm.SetStopGain(mGain);
		cleanAlgorithm.SetAllowNegativeComponents(allowNegative);
			
		std::unique_ptr<AreaSet> cleanAreas;
		if(!cleanAreasFilename.empty())
		{
			cleanAreas.reset(new AreaSet());
			AreaParser parser;
			std::ifstream caFile(cleanAreasFilename.c_str());
			parser.Parse(*cleanAreas, caFile);
			cleanAreas->SetImageProperties(pixelScale, pixelScale, ra, dec, imgWidth, imgHeight);
			cleanAlgorithm.SetCleanAreas(*cleanAreas);
		}
		
		// Start major cleaning loop
		size_t majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			std::cout << " == Cleaning (" << majorIterationNr << ") ==\n";
			cleanAlgorithm.ExecuteMajorIteration(&residual[0], &modelImage[0], &psf[0], imgWidth, imgHeight, reachedMajorThreshold);
			
			if(majorIterationNr == 1)
			{
				std::cout << "Writing residual image... " << std::flush;
				FitsWriter resWriter(std::string(prefixName) + "-residual.fits");
				resWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
				std::cout << "DONE\n";
			}
			
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				FitsWriter modelWriter(std::string(prefixName) + "-model.fits");
				modelWriter.Write(&modelImage[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
				std::cout << "DONE\n";
			}
			
			if(mGain != 1.0)
			{
				std::cout << " == Converting model image to visibilities ==\n";
				inversionAlgorithm->SetAddToModel(false);
				inversionAlgorithm->InvertToVisibilities(&modelImage[0]);
				
				std::cout << " == Constructing image ==\n";
				inversionAlgorithm->SetDoSubtractModel(true);
				inversionAlgorithm->Invert();
				
				if(!reachedMajorThreshold)
				{
					// This was the final major iteration: clean up & save results
					memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
					inversionAlgorithm.reset();
					
					std::cout << "Writing residual image... " << std::flush;
					FitsWriter resWriter(std::string(prefixName) + "-residmajor.fits");
					resWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
					std::cout << "DONE\n";
				}
			}
			
			++majorIterationNr;
			
		} while(reachedMajorThreshold);
		
		std::cout << majorIterationNr << " major iterations were performed.\n";
	}
	
	Model model;
	if(!addModelFilename.empty())
	{
		std::cout << "Reading model from " << addModelFilename << "... " << std::flush;
		model = Model(addModelFilename.c_str());
	}
	CleanAlgorithm::GetModelFromImage(model, &modelImage[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, 0.0, (freqHigh+freqLow)*0.5);
	if(!saveModelFilename.empty())
	{
		std::cout << "Saving model to " << saveModelFilename << "... " << std::flush;
		model.Save(saveModelFilename.c_str());
	}
	
	std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
	ModelRenderer renderer(ra, dec, pixelScale, pixelScale);
	size_t polarizationIndex;
	switch(polarization)
	{
		case InversionAlgorithm::StokesI:
		case InversionAlgorithm::XX: polarizationIndex = 0; break;
		case InversionAlgorithm::XY: polarizationIndex = 1; break;
		case InversionAlgorithm::YX: polarizationIndex = 2; break;
		case InversionAlgorithm::YY: polarizationIndex = 3; break;
	}
	renderer.Restore(&residual[0], imgWidth, imgHeight, model, beamSize, freqLow, freqHigh, polarizationIndex);
	std::cout << "DONE\n";
	
	std::cout << "Writing restored image... " << std::flush;
	FitsWriter restoredWriter(std::string(prefixName) + "-image.fits");
	restoredWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
	std::cout << "DONE\n";
	
}
