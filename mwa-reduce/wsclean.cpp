#include "areaset.h"
#include "beamevaluator.h"
#include "cleanalgorithm.h"
#include "fitswriter.h"
#include "imagebufferallocator.h"
#include "imageweights.h"
#include "inversionalgorithm.h"
#include "modelrenderer.h"
#include "model.h"
#include "msselection.h"
#include "partitionedms.h"
#include "stopwatch.h"
#include "wsinversion.h"

#include "parser/areaparser.h"

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <memory>

std::string commandLine;

void initFitsWriter(FitsWriter& writer, const InversionAlgorithm& inversionAlgorithm, double beamSizeArcmin)
{
	double
		ra = inversionAlgorithm.PhaseCentreRA(),
		dec = inversionAlgorithm.PhaseCentreDec(),
		pixelScaleX = inversionAlgorithm.PixelSizeX(),
		pixelScaleY = inversionAlgorithm.PixelSizeY(),
		freqHigh = inversionAlgorithm.HighestFrequencyChannel(),
		freqLow = inversionAlgorithm.LowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = inversionAlgorithm.BandEnd() - inversionAlgorithm.BandStart(),
		beamSize = inversionAlgorithm.BeamSize(),
		dateObs = inversionAlgorithm.StartTime();
		
	writer.SetImageDimensions(inversionAlgorithm.ImageWidth(), inversionAlgorithm.ImageHeight(), ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(inversionAlgorithm.Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(beamSizeArcmin != 0.0) {
		double beamSizeRad = beamSizeArcmin * (M_PI / 60.0 / 180.0);
		writer.SetBeamInfo(beamSizeRad, beamSizeRad, 0.0);
	}
	else {
		writer.SetBeamInfo(beamSize);
	}
	if(inversionAlgorithm.HasDenormalPhaseCentre())
		writer.SetPhaseCentreShift(inversionAlgorithm.PhaseCentreDL(), inversionAlgorithm.PhaseCentreDM());
	
	writer.SetExtraKeyword("WSCNWLAY", inversionAlgorithm.WGridSize());
	writer.SetExtraKeyword("WSCDATAC", inversionAlgorithm.DataColumnName());
	writer.SetExtraKeyword("WSCWEIGH", inversionAlgorithm.Weighting().ToString());
	writer.SetExtraKeyword("WSCGKRNL", inversionAlgorithm.AntialiasingKernelSize());
	if(inversionAlgorithm.Selection().HasChannelRange())
	{
		writer.SetExtraKeyword("WSCCHANS", inversionAlgorithm.Selection().ChannelRangeStart());
		writer.SetExtraKeyword("WSCCHANE", inversionAlgorithm.Selection().ChannelRangeEnd());
	}
	if(inversionAlgorithm.Selection().HasInterval())
	{
		writer.SetExtraKeyword("WSCTIMES", inversionAlgorithm.Selection().IntervalStart());
		writer.SetExtraKeyword("WSCTIMEE", inversionAlgorithm.Selection().IntervalEnd());
	}
	writer.SetExtraKeyword("WSCFIELD", inversionAlgorithm.Selection().FieldId());
}

void SetCleanParameters(FitsWriter& writer, const CleanAlgorithm& clean)
{
	writer.SetExtraKeyword("WSCNITER", clean.MaxNIter());
	writer.SetExtraKeyword("WSCTHRES", clean.Threshold());
	writer.SetExtraKeyword("WSCGAIN", clean.SubtractionGain());
	writer.SetExtraKeyword("WSCMGAIN", clean.StopGain());
	writer.SetExtraKeyword("WSCNEGCM", clean.AllowNegativeComponents());
	writer.SetExtraKeyword("WSCNEGST", clean.StopOnNegativeComponents());
	writer.SetExtraKeyword("WSCSMPSF", clean.ResizePSF());
}

void UpdateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

int main(int argc, char *argv[])
{
	std::cout << "\n"
		" ** This software package is not public. Please do not use it or distribute it **\n"
		" ** without explicit permission from the author (offringa@gmail.com).          **\n"
		" ** The intension is to make the code public at a later time.                  **\n\n";
	
	if(argc < 2)
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
			"\t-channelsout <count>\n"
			"\t   Splits the bandwidth and makes count nr. of images. Default: 1.\n"
			"\t-field <fieldid>\n"
			"\t   Image the given field id. Default: first field (id 0).\n"
			"\t-weight <weightmode>\n"
			"\t   Weightmode can be: natural, mwa, uniform, briggs. Default: uniform. When using Briggs' weighting,\n"
			"\t   add the robustness parameter, like: \"-weight briggs 0.5\".\n"
			"\t-superweight <factor>\n"
			"\t   Increase the weight gridding box size, similar to Casa's superuniform weighting scheme. Default: 1.0\n"
			"\t   The factor can be rational and can be less than one for subpixel weighting.\n"
			"\t-beamsize <arcmin>\n"
			"\t   Set the FWHM beam size in arcmin for restoring the clean components. Default: longest projected\n"
			"\t   baseline defines restoring beam.\n"
			"\t-makepsf\n"
			"\t   Always make the psf, even when no cleaning is performed.\n"
			"\t-imaginarypart\n"
			"\t   saves the imaginary part instead of the real part; only sensible for xy/yx. Not the default.\n"
			"\t-datacolumn <columnname>\n"
			"\t   Default: CORRECTED_DATA if it exists, otherwise DATA will be used.\n"
			"\t-gkernelsize <size>\n"
			"\t   Gridding antialiasing kernel size. Default: 7.\n"
			"\t-oversampling <factor>\n"
			"\t   Oversampling factor used during gridding. Default: 63.\n"
			"\t-reorder\n"
			"\t-no-reorder\n"
			"\t   Force or disable reordering of Measurement Set. This can be faster when the measurement set needs to\n"
			"\t   be iterated several times, such as with many major iterations or in channel imaging mode.\n"
			"\t   Default: only reorder when in channel imaging mode.\n"
			"\t-addmodel <modelfile>\n"
			"\t-addmodelapp <modelfile>\n"
			"\t-savemodel <modelfile>\n"
			"\t-wlimit <percentage>\n"
			"\t   Do not grid visibilities with a w-value higher than the given percentage of the max w, to save speed\n"
			"\t   (see Tasse et al., 2013, App C),.  Default: grid everything\n"
			"\t-mem <percentage>\n"
			"\t   Limit memory usage to the given fraction of the total system memory. This is an approximate value.\n"
			"\t   Default: 1.\n";
		return -1;
	}
	
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048, channelsOut = 1;
	double pixelScale = 0.01 * M_PI / 180.0, threshold = 0.0, gain = 0.1, mGain = 1.0, beamSize = 0.0, memFraction = 1.0, wLimit = 0.0;
	size_t nWLayers = 0, nIter = 0, antialiasingKernelSize = 7, overSamplingFactor = 63;
	MSSelection globalSelection;
	std::string columnName, addModelFilename, saveModelFilename, cleanAreasFilename;
	PolarizationEnum polarization = Polarization::StokesI;
	WeightMode weightMode(WeightMode::UniformWeighted);
	std::string prefixName = "wsclean";
	bool allowNegative = true, smallPSF = false, addApparentModel = false, stopOnNegative = false, imaginaryPart = false, makePsf = false;
	bool forceReorder = false, forceNoReorder = false;
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
			globalSelection.SetInterval(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "channelrange")
		{
			globalSelection.SetChannelRange(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "channelsout")
		{
			++argi;
			channelsOut = atoi(argv[argi]);
		}
		else if(param == "field")
		{
			++argi;
			globalSelection.SetFieldId(atoi(argv[argi]));
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
		else if(param == "beamsize")
		{
			++argi;
			beamSize = atof(argv[argi]);
		}
		else if(param == "gkernelsize")
		{
			++argi;
			antialiasingKernelSize = atoi(argv[argi]);
		}
		else if(param == "oversampling")
		{
			++argi;
			overSamplingFactor = atoi(argv[argi]);
		}
		else if(param == "reorder")
		{
			forceReorder = true;
			forceNoReorder = false;
		}
		else if(param == "no-reorder")
		{
			forceNoReorder = true;
			forceReorder = false;
		}
		else if(param == "mem")
		{
			++argi;
			memFraction = atof(argv[argi]) / 100.0;
		}
		else if(param == "wlimit")
		{
			++argi;
			wLimit = atof(argv[argi]);
		}
		else {
			throw std::runtime_error("Unknown parameter: " + param);
		}
		
		++argi;
	}
	
	if(argi == argc)
		throw std::runtime_error("No input measurement sets given.");
	
	
	ImageBufferAllocator<double> imageAllocator;
	
	// If no column specified, determine column to use
	if(columnName.empty())
	{
		casa::MeasurementSet ms(argv[argi]);
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			columnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			columnName= "DATA";
		}
	}

	bool doReorder = ((channelsOut != 1) || forceReorder) && !forceNoReorder;
	
	std::vector<PartitionedMS::Handle> partitionedMSHandles;
	if(doReorder)
	{
		for(int i=argi; i != argc; ++i)
		{
			partitionedMSHandles.push_back(PartitionedMS::Partition(argv[i], channelsOut, globalSelection, columnName, true, mGain != 1.0, polarization));
		}
	}
	
	for(size_t outChannelIndex=0; outChannelIndex!=channelsOut; ++outChannelIndex)
	{
		MSSelection partSelection = globalSelection;
		std::unique_ptr<InversionAlgorithm> inversionAlgorithm(new WSInversion(&imageAllocator, memFraction));
		static_cast<WSInversion&>(*inversionAlgorithm).SetGridMode(gridMode);
		
		std::vector<MSProvider*> msProviders;
		for(int i=argi; i != argc; ++i) {
			MSProvider* msProvider;
			if(doReorder)
				msProvider = new PartitionedMS(partitionedMSHandles[i-argi], outChannelIndex);
			else
				msProvider = new ContiguousMS(argv[i], columnName, partSelection, polarization, mGain != 1.0);
			inversionAlgorithm->AddMeasurementSet(msProvider);
			msProviders.push_back(msProvider);
			filenames.push_back(argv[i]);
		}
		if(doReorder)
		{
			size_t startCh, endCh;
			if(globalSelection.HasChannelRange())
			{
				startCh = globalSelection.ChannelRangeStart();
				endCh = globalSelection.ChannelRangeEnd();
			}
			else {
				BandData band(inversionAlgorithm->MeasurementSet(0).MS().spectralWindow());
				startCh = 0;
				endCh = band.ChannelCount();
			}
			size_t newStart = startCh + (endCh - startCh) * outChannelIndex / channelsOut;
			size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / channelsOut;
			partSelection.SetChannelRange(newStart, newEnd);
		}
		std::string partPrefixName;
		if(channelsOut != 1)
		{
			std::ostringstream partPrefixNameStr;
			partPrefixNameStr << prefixName << '-';
			if(outChannelIndex < 1000) partPrefixNameStr << '0';
			if(outChannelIndex < 100) partPrefixNameStr << '0';
			if(outChannelIndex < 10) partPrefixNameStr << '0';
			partPrefixNameStr << outChannelIndex;
			partPrefixName = partPrefixNameStr.str();
		}
		else {
			partPrefixName = prefixName;
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
				imageWeights->Grid(inversionAlgorithm->MeasurementSet(i), weightMode, partSelection);
				if(inversionAlgorithm->MeasurementSetCount() > 1)
					std::cout << i << ' ' << std::flush;
			}
			std::cout << "DONE\n";
			inversionAlgorithm->SetPrecalculatedWeightInfo(imageWeights.get());
		}
		
		inversionAlgorithm->SetImageWidth(imgWidth);
		inversionAlgorithm->SetImageHeight(imgHeight);
		inversionAlgorithm->SetPixelSizeX(pixelScale);
		inversionAlgorithm->SetPixelSizeY(pixelScale);
		if(nWLayers != 0)
			inversionAlgorithm->SetWGridSize(nWLayers);
		else
			inversionAlgorithm->SetNoWGridSize();
		inversionAlgorithm->SetAntialiasingKernelSize(antialiasingKernelSize);
		inversionAlgorithm->SetOverSamplingFactor(overSamplingFactor);
		inversionAlgorithm->SetPolarization(polarization);
		inversionAlgorithm->SetDataColumnName(columnName);
		inversionAlgorithm->SetWeighting(weightMode);
		inversionAlgorithm->SetImaginaryPart(imaginaryPart);
		inversionAlgorithm->SetSelection(partSelection);
		inversionAlgorithm->SetWLimit(wLimit/100.0);
		
		double* psf = 0;
		bool isFirstInversion = true;
		
		Stopwatch inversionWatch(false), predictingWatch(false), cleaningWatch(false);
		
		if(nIter > 0 || makePsf)
		{
			std::cout << std::flush << " == Constructing PSF ==\n";
			inversionWatch.Start();
			inversionAlgorithm->SetDoImagePSF(true);
			inversionAlgorithm->SetVerbose(isFirstInversion);
			inversionAlgorithm->Invert();
				
			psf = imageAllocator.Allocate(imgWidth * imgHeight);
			memcpy(psf, inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
			inversionWatch.Pause();
			
			isFirstInversion = false;
			std::cout << "Beam size is " << inversionAlgorithm->BeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
			
			std::cout << "Writing psf image... " << std::flush;
			FitsWriter fitsWriter;
			initFitsWriter(fitsWriter, *inversionAlgorithm, beamSize);
			fitsWriter.Write(std::string(partPrefixName) + "-psf.fits", psf);
			std::cout << "DONE\n";
			
			CleanAlgorithm::RemoveNaNsInPSF(psf, imgWidth, imgHeight);
		}
		
		std::cout << std::flush << " == Constructing image ==\n";
		inversionWatch.Start();
		if(nWLayers != 0)
			inversionAlgorithm->SetWGridSize(nWLayers);
		else
			inversionAlgorithm->SetNoWGridSize();
		inversionAlgorithm->SetDoImagePSF(false);
		inversionAlgorithm->SetVerbose(isFirstInversion);
		inversionAlgorithm->Invert();
		inversionWatch.Pause();
		inversionAlgorithm->SetVerbose(false);
		
		if(inversionAlgorithm->HasGriddingCorrectionImage())
		{
			std::cout << "Writing gridding correction image... " << std::flush;
			double* gridding = imageAllocator.Allocate(imgWidth * imgHeight);
			inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
			FitsWriter fitsWriter;
			initFitsWriter(fitsWriter, *inversionAlgorithm, beamSize);
			fitsWriter.Write(std::string(partPrefixName) + "-gridding.fits", &gridding[0]);
			imageAllocator.Free(gridding);
			std::cout << "DONE\n";
		}
			
		double
			*modelImage = imageAllocator.Allocate(imgWidth * imgHeight),
			*residual = imageAllocator.Allocate(imgWidth * imgHeight);
		memcpy(residual, inversionAlgorithm->ImageResult(), imgWidth * imgHeight * sizeof(double));
		memset(modelImage, 0, imgWidth * imgHeight * sizeof(double));
		
		std::cout << "Writing dirty image... " << std::flush;
		FitsWriter fitsWriter;
		initFitsWriter(fitsWriter, *inversionAlgorithm, beamSize);
		fitsWriter.Write(std::string(partPrefixName) + "-dirty.fits", residual);
		std::cout << "DONE\n";
		
		if(mGain == 1.0)
			inversionAlgorithm.reset();
		
		CleanAlgorithm cleanAlgorithm;
		cleanAlgorithm.SetMaxNIter(nIter);
		cleanAlgorithm.SetThreshold(threshold);
		cleanAlgorithm.SetSubtractionGain(gain);
		cleanAlgorithm.SetStopGain(mGain);
		cleanAlgorithm.SetAllowNegativeComponents(allowNegative);
		cleanAlgorithm.SetStopOnNegativeComponents(stopOnNegative);
		cleanAlgorithm.SetResizePSF(smallPSF);

		SetCleanParameters(fitsWriter, cleanAlgorithm);
		UpdateCleanParameters(fitsWriter, 0, 0);
			
		if(nIter > 0)
		{
			std::unique_ptr<AreaSet> cleanAreas;
			if(!cleanAreasFilename.empty())
			{
				cleanAreas.reset(new AreaSet());
				AreaParser parser;
				std::ifstream caFile(cleanAreasFilename.c_str());
				parser.Parse(*cleanAreas, caFile);
				cleanAreas->SetImageProperties(pixelScale, pixelScale, inversionAlgorithm->PhaseCentreRA(), inversionAlgorithm->PhaseCentreDec(), imgWidth, imgHeight);
				cleanAlgorithm.SetCleanAreas(*cleanAreas);
			}
			
			// Start major cleaning loop
			size_t majorIterationNr = 1;
			bool reachedMajorThreshold = false;
			do {
				std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
				cleaningWatch.Start();
				cleanAlgorithm.ExecuteMajorIteration(residual, modelImage, psf, imgWidth, imgHeight, reachedMajorThreshold);
				cleaningWatch.Pause();
				
				UpdateCleanParameters(fitsWriter, cleanAlgorithm.IterationNumber(), majorIterationNr);
				
				if(majorIterationNr == 1)
				{
					if(mGain == 1.0)
					{
						std::cout << "Writing residual image... " << std::flush;
						fitsWriter.Write(std::string(partPrefixName) + "-residual.fits", &residual[0]);
					}
					else {
						std::cout << "Writing first iteration image... " << std::flush;
						fitsWriter.Write(std::string(partPrefixName) + "-first-residual.fits", &residual[0]);
					}
					std::cout << "DONE\n";
				}
				
				if(!reachedMajorThreshold)
				{
					std::cout << "Writing model image... " << std::flush;
					fitsWriter.Write(std::string(partPrefixName) + "-model.fits", &modelImage[0]);
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
						fitsWriter.Write(std::string(partPrefixName) + "-residual.fits", &residual[0]);
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
			default: throw std::runtime_error("Unsupported polarization");
		}
		double
			freqLow = fitsWriter.Frequency() - fitsWriter.Bandwidth()*0.5,
			freqHigh = fitsWriter.Frequency() + fitsWriter.Bandwidth()*0.5;
		renderer.Restore(&residual[0], imgWidth, imgHeight, model, fitsWriter.BeamSizeMajorAxis(), freqLow, freqHigh, polarizationIndex);
		std::cout << "DONE\n";
		
		std::cout << "Writing restored image... " << std::flush;
		fitsWriter.Write(std::string(partPrefixName) + "-image.fits", &residual[0]);
		std::cout << "DONE\n";
		
		imageAllocator.ReportStatistics();
		std::cout << "Inversion: " << inversionWatch.ToString() << ", prediction: " << predictingWatch.ToString() << ", cleaning: " << cleaningWatch.ToString() << '\n';
		
		imageAllocator.Free(residual);
		imageAllocator.Free(modelImage);
		imageAllocator.Free(psf);
		
		for(std::vector<MSProvider*>::iterator i=msProviders.begin(); i != msProviders.end(); ++i)
			delete *i;
	}
}
