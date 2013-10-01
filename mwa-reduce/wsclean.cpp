#include "cleanalgorithm.h"
#include "inversionalgorithm.h"
#include "wsinversion.h"
#include "fitswriter.h"
#include "modelrenderer.h"
#include "model.h"
#include "areaset.h"
#include "parser/areaparser.h"
#include "beamevaluator.h"
#include "imageweights.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>

int main(int argc, char *argv[])
{
	std::cout << "\n"
		" ** This software package is not public. Please do not use it or distribute it **\n"
		" ** without explicit permission from the author (offringa@gmail.com).          **\n"
		" ** The intension is to make the code public at a later time.                  **\n\n";
	
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
			"\t-smallpsf\n"
			"\t   Resize the psf to speed up minor clean iterations.\n"
			"\t-pol <xx, yy, xy, yx or stokesi>\n"
			"\t-negative\n"
			"\t   Allow negative components during cleaning\n"
			"\t-interval <start-index> <stop-index>\n"
			"\t   Only image the given interval. Indices specify the timesteps, stop is exclusive.\n"
			"\t-weight <weightmode>\n"
			"\t   weightmode can be: natural, mwa, uniform\n"
			"\t-datacolumn <columnname>\n"
			"\t-addmodel <modelfile>\n"
			"\t-addmodelapp <modelfile>\n"
			"\t-savemodel <modelfile>\n";
		return -1;
	}
	
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0, threshold = 0.0, gain = 0.1, mGain = 1.0;
	size_t nWLayers = 64, nIter = 0, intervalStart = 0, intervalStop = 0;
	std::string columnName, addModelFilename, saveModelFilename, cleanAreasFilename;
	enum InversionAlgorithm::PolarizationEnum polarization = InversionAlgorithm::StokesI;
	enum InversionAlgorithm::WeightingEnum weightMode = InversionAlgorithm::DistanceWeighted;
	std::string prefixName = "wsclean";
	bool allowNegative = false, smallPSF = false, addApparentModel = false;
	enum LayeredImager::GridModeEnum gridMode = LayeredImager::NearestNeighbour;
	std::vector<std::string> filenames;
	
	while(argi < argc && argv[argi][0] == '-')
	{
		const std::string param = &argv[argi][1];
		if(param == "size")
		{
			imgWidth = atoi(argv[argi+1]);
			imgHeight = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(param == "scale")
		{
			++argi;
			pixelScale = atof(argv[argi]) * M_PI / 180.0;
		}
		else if(param == "nwlayers")
		{
			++argi;
			nWLayers = atoi(argv[argi]);
		}
		else if(param == "gain")
		{
			++argi;
			gain = atof(argv[argi]);
		}
		else if(param == "mgain")
		{
			++argi;
			mGain = atof(argv[argi]);
		}
		else if(param == "niter")
		{
			++argi;
			nIter = atoi(argv[argi]);
		}
		else if(param == "threshold")
		{
			++argi;
			threshold = atof(argv[argi]);
		}
		else if(param == "datacolumn")
		{
			++argi;
			columnName = argv[argi];
		}
		else if(param == "pol")
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
		else if(param == "negative")
		{
			allowNegative = true;
		}
		else if(param == "addmodel")
		{
			++argi;
			addModelFilename = argv[argi];
		}
		else if(param == "addmodelapp")
		{
			++argi;
			addModelFilename = argv[argi];
			addApparentModel = true;
		}
		else if(param == "savemodel")
		{
			++argi;
			saveModelFilename = argv[argi];
		}
		else if(param == "cleanareas")
		{
			++argi;
			cleanAreasFilename = argv[argi];
		}
		else if(param == "name")
		{
			++argi;
			prefixName = argv[argi];
		}
		else if(param == "gridmode")
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
		else if(param == "smallpsf")
		{
			smallPSF = true;
		}
		else if(param == "interval")
		{
			intervalStart = atoi(argv[argi+1]);
			intervalStop = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(param == "weight")
		{
			++argi;
			std::string weightArg = argv[argi];
			if(weightArg == "natural")
				weightMode = InversionAlgorithm::NaturalWeighted;
			else if(weightArg == "mwa")
				weightMode = InversionAlgorithm::DistanceWeighted;
			else if(weightArg == "uniform")
				weightMode = InversionAlgorithm::UniformWeighted;
		}
		else {
			throw std::runtime_error("Unknown parameter: " + param);
		}
		
		++argi;
	}
	
	if(argi == argc)
		throw std::runtime_error("No input measurement sets given.");
	
	std::unique_ptr<InversionAlgorithm> inversionAlgorithm(new WSInversion());
	static_cast<WSInversion&>(*inversionAlgorithm).SetGridMode(gridMode);
	
	for(int i=argi; i != argc; ++i) {
		inversionAlgorithm->AddMeasurementSetPath(argv[i]);
		filenames.push_back(argv[i]);
	}
	
	// If uniform weighted; initialize weight grid.
	std::unique_ptr<ImageWeights> imageWeights;
	if(weightMode == InversionAlgorithm::UniformWeighted)
	{
		std::cout << "Precalculating weights for uniform weighting... " << std::flush;
		imageWeights.reset(new ImageWeights(imgWidth, imgHeight, pixelScale));
		for(size_t i=0; i!=inversionAlgorithm->MeasurementSetCount(); ++i)
		{
			casa::MeasurementSet ms(inversionAlgorithm->MeasurementSetPath(i));
			imageWeights->Grid(ms);
			if(inversionAlgorithm->MeasurementSetCount() > 1)
				std::cout << i << ' ' << std::flush;
		}
		std::cout << "DONE\n";
		inversionAlgorithm->SetPrecalculatedWeightInfo(imageWeights.get());
	}
	
	// If no column specified, determine column to use
	if(columnName.empty())
	{
		casa::MeasurementSet ms(inversionAlgorithm->MeasurementSetPath(0));
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		std::string dataColumn;
		if(hasCorrected) {
			std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			columnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			columnName= "DATA";
		}
	}
	
	inversionAlgorithm->SetImageWidth(imgWidth);
	inversionAlgorithm->SetImageHeight(imgHeight);
	inversionAlgorithm->SetPixelSizeX(pixelScale);
	inversionAlgorithm->SetPixelSizeY(pixelScale);
	inversionAlgorithm->SetWGridSize(nWLayers);
	inversionAlgorithm->SetPolarization(polarization);
	inversionAlgorithm->SetDataColumnName(columnName);
	inversionAlgorithm->SetWeighting(weightMode);
	if(intervalStop != 0)
	{
		inversionAlgorithm->SetInterval(intervalStart, intervalStop);
	}
	
	double
		ra, dec,
		freqHigh, freqLow, freqCentre,
		bandwidth, beamSize, dateObs;
		
	std::vector<double> psf;
	
	if(nIter > 0)
	{
		std::cout << std::flush << " == Constructing PSF ==\n";
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
		std::cout << "Beam size is " << beamSize*(180.0*60.0/M_PI) << " arcmin.\n";
		
		std::cout << "Writing psf image... " << std::flush;
		FitsWriter psfWriter(std::string(prefixName) + "-psf.fits");
		psfWriter.SetBeamInfo(beamSize);
		psfWriter.Write(&psf[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
		std::cout << "DONE\n";
		
		CleanAlgorithm::RemoveNaNsInPSF(&psf[0], imgWidth, imgHeight);
	
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
	
	std::cout << std::flush << " == Constructing image ==\n";
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
	dirtyWriter.SetBeamInfo(beamSize);
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
		cleanAlgorithm.SetResizePSF(smallPSF);
			
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
			std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
			cleanAlgorithm.ExecuteMajorIteration(&residual[0], &modelImage[0], &psf[0], imgWidth, imgHeight, reachedMajorThreshold);
			
			if(majorIterationNr == 1)
			{
				std::cout << "Writing residual image... " << std::flush;
				FitsWriter resWriter(std::string(prefixName) + "-residual.fits");
				resWriter.SetBeamInfo(beamSize);
				resWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
				std::cout << "DONE\n";
			}
			
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				FitsWriter modelWriter(std::string(prefixName) + "-model.fits");
				modelWriter.SetBeamInfo(beamSize);
				modelWriter.Write(&modelImage[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
				std::cout << "DONE\n";
			}
			
			if(mGain != 1.0)
			{
				std::cout << std::flush << " == Converting model image to visibilities ==\n";
				inversionAlgorithm->SetAddToModel(false);
				inversionAlgorithm->InvertToVisibilities(&modelImage[0]);
				
				std::cout << std::flush << " == Constructing image ==\n";
				inversionAlgorithm->SetDoSubtractModel(true);
				inversionAlgorithm->Invert();
				
				memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
				
				if(!reachedMajorThreshold)
				{
					// This was the final major iteration: clean up & save results
					inversionAlgorithm.reset();
					
					std::cout << "Writing residual image... " << std::flush;
					FitsWriter resWriter(std::string(prefixName) + "-residmajor.fits");
					resWriter.SetBeamInfo(beamSize);
					resWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
					std::cout << "DONE\n";
				}
				
				++majorIterationNr;
			}
			
		} while(reachedMajorThreshold);
		
		std::cout << majorIterationNr << " major iterations were performed.\n";
	}
	
	Model model;
	if(!addModelFilename.empty())
	{
		std::cout << "Reading model from " << addModelFilename << "... " << std::flush;
		model = Model(addModelFilename.c_str());
		if(addApparentModel)
		{
			casa::MeasurementSet ms(filenames[0]);
			BeamEvaluator beamEval(ms, false);
			beamEval.AbsToApparent(model);
		}
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
	restoredWriter.SetBeamInfo(beamSize);
	restoredWriter.Write(&residual[0], imgWidth, imgHeight, ra, dec, pixelScale, pixelScale, freqCentre, bandwidth, dateObs);
	std::cout << "DONE\n";
	
}
