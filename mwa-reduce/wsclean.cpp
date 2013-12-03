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
#include "stopwatch.h"

#include <iostream>
#include <memory>

#include <boost/algorithm/string.hpp>

std::string commandLine;

void initFitsWriter(FitsWriter& writer, const InversionAlgorithm& inversionAlgorithm)
{
	double
		ra = inversionAlgorithm.ImageResultRA(),
		dec = inversionAlgorithm.ImageResultDec(),
		pixelScaleX = inversionAlgorithm.PixelSizeX(),
		pixelScaleY = inversionAlgorithm.PixelSizeY(),
		freqHigh = inversionAlgorithm.ImageHighestFrequencyChannel(),
		freqLow = inversionAlgorithm.ImageLowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = inversionAlgorithm.ImageBandEnd() - inversionAlgorithm.ImageBandStart(),
		beamSize = inversionAlgorithm.ImageBeamSize(),
		dateObs = inversionAlgorithm.ImageStartTime();
		
	writer.SetImageDimensions(inversionAlgorithm.ImageWidth(), inversionAlgorithm.ImageHeight(), ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetBeamInfo(beamSize);
	writer.SetPolarization(inversionAlgorithm.Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
}

int main(int argc, char *argv[])
{
	std::cout << "\n"
		" ** This software package is not public. Please do not use it or distribute it **\n"
		" ** without explicit permission from the author (offringa@gmail.com).          **\n"
		" ** The intension is to make the code public at a later time.                  **\n\n";
	
	if(argc < 3)
	{
		std::cout << "Syntax:\twsclean [options] <input-ms> [<2nd-ms> [..]]\n"
			"Will create cleaned images of the input ms(es).\n"
			"If multiple mses are specified, they need to be phase-rotated to the same point on the sky.\n"
			"Options can be:\n"
			"\t-name <image-prefix>\n"
			"\t   Use image-prefix as prefix for output files. Default is 'wsclean'.\n"
			"\t-size <width> <height>\n"
			"\t   Default: 2048 x 2048\n"
			"\t-scale <pixel-scale>\n"
			"\t   Scale of a pixel in degrees, e.g. 0.012. Default: 0.01\n"
			"\t-nwlayers <nwlayers>\n"
			"\t   Number of w-layers to use. Default: minimum suggested #w-layers for first MS.\n"
			"\t-niter <niter>\n"
			"\t   Maximum number of clean iterations to perform. Default: 0\n"
			"\t-threshold <threshold>\n"
			"\t   Stopping clean thresholding in Jy. Default: 0.0\n"
			"\t-gain <gain>\n"
			"\t   Cleaning gain: Ratio of peak that will be subtracted in each iteration. Default: 0.1\n"
			"\t-mgain <gain>\n"
			"\t   Cleaning gain for major iterations: Ratio of peak that will be subtracted in each major\n"
			"\t   iteration (default = 1.0, to use major iterations, 0.9 is a good value). Default: 1.0\n"
			"\t-smallpsf\n"
			"\t   Resize the psf to speed up minor clean iterations. Not the default.\n"
			"\t-pol <xx, yy, xy, yx or stokesi>\n"
			"\t   Default: stokesi.\n"
			"\t-gridmode <nn or kb>\n"
			"\t   Kernel and mode used for gridding: kb = Kaiser-Bessel (currently 7 pixels), nn = nearest\n"
			"\t   neighbour (no kernel). Default: kb.\n"
			"\t-nonegative\n"
			"\t   Do not allow negative components during cleaning. Not the default.\n"
			"\t-negative\n"
			"\t   Default on: opposite of -nonegative.\n"
			"\t-stopnegative\n"
			"\t   Stop on negative components. Not the default.\n"
			"\t-interval <start-index> <end-index>\n"
			"\t   Only image the given time interval. Indices specify the timesteps, end index is exclusive.\n"
			"\t   Default: image all time steps.\n"
			"\t-channelrange <start-channel> <end-channel>\n"
			"\t   Only image the given channel range. Indices specify channel indices, end index is exclusive.\n"
			"\t   Default: image all channels.\n"
			"\t-weight <weightmode>\n"
			"\t   Weightmode can be: natural, mwa, uniform, briggs. Default: uniform. When using Briggs' weighting,\n"
			"\t   add the robustness parameter, like: \"-weight briggs 0.5\".\n"
			"\t-superweight <factor>\n"
			"\t   Increase the weight gridding box size, similar to Casa's superuniform weighting scheme. Default: 1.0\n"
			"\t   The factor can be rational and can be less than one for subpixel weighting.\n"
			"\t-makepsf\n"
			"\t   Always make the psf, even when no cleaning is performed.\n"
			"\t-imaginarypart\n"
			"\t   saves the imaginary part instead of the real part; only sensible for xy/yx. Not the default.\n"
			"\t-datacolumn <columnname>\n"
			"\t   Default: CORRECTED_DATA if it exists, otherwise DATA will be used.\n"
			"\t-addmodel <modelfile>\n"
			"\t-addmodelapp <modelfile>\n"
			"\t-savemodel <modelfile>\n";
		return -1;
	}
	
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0, threshold = 0.0, gain = 0.1, mGain = 1.0;
	size_t nWLayers = 0, nIter = 0, intervalStart = 0, intervalEnd = 0, channelRangeStart = 0, channelRangeEnd = 0;
	std::string columnName, addModelFilename, saveModelFilename, cleanAreasFilename;
	PolarizationEnum polarization = Polarization::StokesI;
	WeightMode weightMode(WeightMode::UniformWeighted);
	std::string prefixName = "wsclean";
	bool allowNegative = true, smallPSF = false, addApparentModel = false, stopOnNegative = false, imaginaryPart = false, makePsf = false;
	enum LayeredImager::GridModeEnum gridMode = LayeredImager::KaiserBessel;
	std::vector<std::string> filenames;
	
	while(argi < argc && argv[argi][0] == '-')
	{
		const std::string param = &argv[argi][1];
		if(param == "size")
		{
			imgWidth = atoi(argv[argi+1]);
			imgHeight = atoi(argv[argi+2]);
			if(imgWidth != imgHeight)
				throw std::runtime_error("width != height : Can't handle non-square images yet");
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
				polarization = Polarization::XX;
			else if(polStr == "xy")
				polarization = Polarization::XY;
			else if(polStr == "yx")
				polarization = Polarization::YX;
			else if(polStr == "yy")
				polarization = Polarization::YY;
			else if(polStr == "stokesi" || polStr == "i")
				polarization = Polarization::StokesI;
			else
				throw std::runtime_error("Unknown polarization given in -pol parameter");
		}
		else if(param == "imaginarypart")
		{
			imaginaryPart = true;
		}
		else if(param == "negative")
		{
			allowNegative = true;
		}
		else if(param == "nonegative")
		{
			allowNegative = false;
		}
		else if(param == "stopnegative")
		{
			stopOnNegative = true;
		}
		else if(param == "makepsf")
		{
			makePsf = true;
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
			intervalEnd = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(param == "channelrange")
		{
			channelRangeStart = atoi(argv[argi+1]);
			channelRangeEnd = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(param == "weight")
		{
			++argi;
			std::string weightArg = argv[argi];
			if(weightArg == "natural")
				weightMode.SetMode(WeightMode(WeightMode::NaturalWeighted));
			else if(weightArg == "mwa")
				weightMode.SetMode(WeightMode(WeightMode::DistanceWeighted));
			else if(weightArg == "uniform")
				weightMode.SetMode(WeightMode(WeightMode::UniformWeighted));
			else if(weightArg == "briggs")
			{
				++argi;
				double robustness = atof(argv[argi]);
				weightMode.SetMode(WeightMode::Briggs(robustness));
			}
			else throw std::runtime_error("Unknown weighting mode specified");
		}
		else if(param == "superweight")
		{
			++argi;
			weightMode.SetSuperWeight(atof(argv[argi]));
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
	
	// Store command line, to write later in fits file
	std::ostringstream commandLineStr;
	commandLineStr << "wsclean";
	for(int i=1; i!=argc; ++i)
		commandLineStr << ' ' << argv[i];
	commandLine = commandLineStr.str();
	
	// Initialize weight grid if necessary.
	std::unique_ptr<ImageWeights> imageWeights;
	if(weightMode.RequiresGridding())
	{
		std::cout << "Precalculating weights for " << weightMode.ToString() << " weighting... " << std::flush;
		imageWeights.reset(new ImageWeights(imgWidth, imgHeight, pixelScale, weightMode.SuperWeight()));
		for(size_t i=0; i!=inversionAlgorithm->MeasurementSetCount(); ++i)
		{
			casa::MeasurementSet ms(inversionAlgorithm->MeasurementSetPath(i));
			imageWeights->Grid(ms, weightMode);
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
	if(nWLayers != 0)
		inversionAlgorithm->SetWGridSize(nWLayers);
	else
		inversionAlgorithm->SetNoWGridSize();
	inversionAlgorithm->SetPolarization(polarization);
	inversionAlgorithm->SetDataColumnName(columnName);
	inversionAlgorithm->SetWeighting(weightMode);
	inversionAlgorithm->SetImaginaryPart(imaginaryPart);
	if(intervalEnd != 0)
		inversionAlgorithm->SetInterval(intervalStart, intervalEnd);
	if(channelRangeEnd != 0)
		inversionAlgorithm->SetChannelRange(channelRangeStart, channelRangeEnd);
	
	std::vector<double> psf;
	bool isFirstInversion = true;
	
	Stopwatch inversionWatch(false), predictingWatch(false), cleaningWatch(false);
	
	if(nIter > 0 || makePsf)
	{
		std::cout << std::flush << " == Constructing PSF ==\n";
		inversionWatch.Start();
		inversionAlgorithm->SetDoImagePSF(true);
		inversionAlgorithm->SetVerbose(isFirstInversion);
		inversionAlgorithm->Invert();
		psf.resize(imgWidth * imgHeight);
		memcpy(&psf[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
		inversionWatch.Pause();
		
		isFirstInversion = false;
		std::cout << "Beam size is " << inversionAlgorithm->ImageBeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
		
		std::cout << "Writing psf image... " << std::flush;
		FitsWriter fitsWriter;
		initFitsWriter(fitsWriter, *inversionAlgorithm);
		fitsWriter.Write(std::string(prefixName) + "-psf.fits", &psf[0]);
		std::cout << "DONE\n";
		
		CleanAlgorithm::RemoveNaNsInPSF(&psf[0], imgWidth, imgHeight);
	
		if(inversionAlgorithm->HasGriddingCorrectionImage())
		{
			std::cout << "Writing gridding correction image... " << std::flush;
			std::vector<double> gridding(imgWidth * imgHeight);
			inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
			fitsWriter.Write(std::string(prefixName) + "-gridding.fits", &gridding[0]);
			std::cout << "DONE\n";
		}
	}
	
	std::cout << std::flush << " == Constructing image ==\n";
	inversionWatch.Start();
	inversionAlgorithm->SetDoImagePSF(false);
	inversionAlgorithm->SetVerbose(isFirstInversion);
	inversionAlgorithm->Invert();
	inversionWatch.Pause();
	inversionAlgorithm->SetVerbose(false);
	
	std::vector<double> modelImage(imgWidth * imgHeight), residual(imgWidth * imgHeight);
	memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
	
	std::cout << "Writing dirty image... " << std::flush;
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter, *inversionAlgorithm);
	fitsWriter.Write(std::string(prefixName) + "-dirty.fits", &residual[0]);
	std::cout << "DONE\n";
	
	if(mGain == 1.0)
		inversionAlgorithm.reset();
	
	if(nIter > 0)
	{
		CleanAlgorithm cleanAlgorithm;
		cleanAlgorithm.SetMaxNIter(nIter);
		cleanAlgorithm.SetThreshold(threshold);
		cleanAlgorithm.SetSubtractionGain(gain);
		cleanAlgorithm.SetStopGain(mGain);
		cleanAlgorithm.SetAllowNegativeComponents(allowNegative);
		cleanAlgorithm.SetStopOnNegativeComponents(stopOnNegative);
		cleanAlgorithm.SetResizePSF(smallPSF);
			
		std::unique_ptr<AreaSet> cleanAreas;
		if(!cleanAreasFilename.empty())
		{
			cleanAreas.reset(new AreaSet());
			AreaParser parser;
			std::ifstream caFile(cleanAreasFilename.c_str());
			parser.Parse(*cleanAreas, caFile);
			cleanAreas->SetImageProperties(pixelScale, pixelScale, inversionAlgorithm->ImageResultRA(), inversionAlgorithm->ImageResultDec(), imgWidth, imgHeight);
			cleanAlgorithm.SetCleanAreas(*cleanAreas);
		}
		
		// Start major cleaning loop
		size_t majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
			cleaningWatch.Start();
			cleanAlgorithm.ExecuteMajorIteration(&residual[0], &modelImage[0], &psf[0], imgWidth, imgHeight, reachedMajorThreshold);
			cleaningWatch.Pause();
			
			if(majorIterationNr == 1)
			{
				std::cout << "Writing residual image... " << std::flush;
				fitsWriter.Write(std::string(prefixName) + "-residual.fits", &residual[0]);
				std::cout << "DONE\n";
			}
			
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				fitsWriter.Write(std::string(prefixName) + "-model.fits", &modelImage[0]);
				std::cout << "DONE\n";
			}
			
			if(mGain != 1.0)
			{
				std::cout << std::flush << " == Converting model image to visibilities ==\n";
				predictingWatch.Start();
				inversionAlgorithm->SetAddToModel(false);
				inversionAlgorithm->InvertToVisibilities(&modelImage[0]);
				predictingWatch.Pause();
				
				std::cout << std::flush << " == Constructing image ==\n";
				inversionWatch.Start();
				inversionAlgorithm->SetDoSubtractModel(true);
				inversionAlgorithm->Invert();
				inversionWatch.Pause();
				
				memcpy(&residual[0], inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
				
				if(!reachedMajorThreshold)
				{
					// This was the final major iteration: clean up & save results
					inversionAlgorithm.reset();
					
					std::cout << "Writing residual image... " << std::flush;
					fitsWriter.Write(std::string(prefixName) + "-residmajor.fits", &residual[0]);
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
	CleanAlgorithm::GetModelFromImage(model, &modelImage[0], imgWidth, imgHeight, fitsWriter.RA(), fitsWriter.Dec(), pixelScale, pixelScale, 0.0, fitsWriter.Frequency(), polarization);
	if(!saveModelFilename.empty())
	{
		std::cout << "Saving model to " << saveModelFilename << "... " << std::flush;
		model.Save(saveModelFilename.c_str());
	}
	
	std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
	ModelRenderer renderer(fitsWriter.RA(), fitsWriter.Dec(), pixelScale, pixelScale);
	size_t polarizationIndex;
	switch(polarization)
	{
		case Polarization::StokesI:
		case Polarization::XX: polarizationIndex = 0; break;
		case Polarization::XY: polarizationIndex = 1; break;
		case Polarization::YX: polarizationIndex = 2; break;
		case Polarization::YY: polarizationIndex = 3; break;
	}
	double
		freqLow = fitsWriter.Frequency() - fitsWriter.Bandwidth()*0.5,
		freqHigh = fitsWriter.Frequency() + fitsWriter.Bandwidth()*0.5;
	renderer.Restore(&residual[0], imgWidth, imgHeight, model, fitsWriter.BeamSizeMajorAxis(), freqLow, freqHigh, polarizationIndex);
	std::cout << "DONE\n";
	
	std::cout << "Writing restored image... " << std::flush;
	fitsWriter.Write(std::string(prefixName) + "-image.fits", &residual[0]);
	std::cout << "DONE\n";
	
	std::cout << "Inversion: " << inversionWatch.ToString() << ", prediction: " << predictingWatch.ToString() << ", cleaning: " << cleaningWatch.ToString() << '\n';
}
