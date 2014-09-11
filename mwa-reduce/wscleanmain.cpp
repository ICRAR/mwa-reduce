#include "wsclean.h"
#include "wscversion.h"

#include <boost/algorithm/string.hpp>

#include <iostream>

int main(int argc, char *argv[])
{
	std::cout << "\n"
		"WSClean version " WSCLEAN_VERSION_STR " (" WSCLEAN_VERSION_DATE ")\n"
		"This software package is released under the GPL version 3.\n"
	  "Author: André Offringa (offringa@gmail.com).\n\n";
#ifndef NDEBUG
	std::cout << "\n"
		"WARNING: Symbol NDEBUG was not defined; this WSClean version was\n"
		"compiled as a DEBUG version. This can seriously affect performance!\n\n";
#endif
	
	if(argc < 2)
	{
		std::cout << "Syntax: wsclean [options] <input-ms> [<2nd-ms> [..]]\n"
			"Will create cleaned images of the input ms(es).\n"
			"If multiple mses are specified, they need to be phase-rotated to the same point on the sky.\n\n"
			"Options can be:\n\n"
			"  ** GENERAL OPTIONS **\n"
			"-version\n"
			"   Print WSClean's version and exit.\n"
			"\n"
			"  ** INVERSION OPTIONS **\n"
			"-name <image-prefix>\n"
			"   Use image-prefix as prefix for output files. Default is 'wsclean'.\n"
			"-size <width> <height>\n"
			"   Default: 2048 x 2048\n"
			"-scale <pixel-scale>\n"
			"   Scale of a pixel in degrees, e.g. 0.012. Default: 0.01\n"
			"-nwlayers <nwlayers>\n"
			"   Number of w-layers to use. Default: minimum suggested #w-layers for first MS.\n"
			"-channelsout <count>\n"
			"   Splits the bandwidth and makes count nr. of images. Default: 1.\n"
			"-predict <image-prefix>\n"
			"   Only perform a single prediction for an existing image. Doesn't do any imaging or cleaning.\n"
			"-nosmallinversion and -smallinversion\n"
			"   Perform inversion at the Nyquist resolution and upscale the image to the requested image size afterwards.\n"
			"   This speeds up inversion considerably, but makes aliasing slightly slightly worse. This effect is\n"
			"   in most cases <1%. Default: on.\n"
			"-weight <weightmode>\n"
			"   Weightmode can be: natural, mwa, uniform, briggs. Default: uniform. When using Briggs' weighting,\n"
			"   add the robustness parameter, like: \"-weight briggs 0.5\".\n"
			"-superweight <factor>\n"
			"   Increase the weight gridding box size, similar to Casa's superuniform weighting scheme. Default: 1.0\n"
			"   The factor can be rational and can be less than one for subpixel weighting.\n"
			"-mfsweighting\n"
			"   In spectral mode, calculate the weights as if the image was made using MFS. This makes sure that the sum of\n"
			"   channel images equals the MFS weights. Otherwise, the channel image will become a bit more naturally weighted.\n"
			"   This is only relevant for weighting modes that require gridding (i.e., Uniform, Briggs').\n"
			"   Default: off, unless -joinchannels is specified.\n"
			"-nomfsweighting\n"
			"   Opposite of -mfsweighting; can be used to turn off MFS weighting in -joinchannels mode.\n"
			"-gridmode <nn or kb>\n"
			"   Kernel and mode used for gridding: kb = Kaiser-Bessel (default with 7 pixels), nn = nearest\n"
			"   neighbour (no kernel). Default: kb.\n"
			"-gkernelsize <size>\n"
			"   Gridding antialiasing kernel size. Default: 7.\n"
			"-oversampling <factor>\n"
			"   Oversampling factor used during gridding. Default: 63.\n"
			"-reorder\n"
			"-no-reorder\n"
			"   Force or disable reordering of Measurement Set. This can be faster when the measurement set needs to\n"
			"   be iterated several times, such as with many major iterations or in channel imaging mode.\n"
			"   Default: only reorder when in channel imaging mode.\n"
			"-makepsf\n"
			"   Always make the psf, even when no cleaning is performed.\n"
			"-j <threads>\n"
			"   Specify number of computing threads to use, i.e., number of cpu cores that will be used.\n"
			"   Default: use all cpu cores.\n"
			"-mem <percentage>\n"
			"   Limit memory usage to the given fraction of the total system memory. This is an approximate value.\n"
			"   Default: 100.\n"
			"-absmem <memory limit>\n"
			"   Like -mem, but this specifies a fixed amount of memory in gigabytes.\n"
			"\n"
			"  ** DATA SELECTION OPTIONS **\n"
			"-pol <list>\n"
			"   Default: \'I\'. Possible values: XX, XY, YX, YY, I, Q, U, V, RR, RL, LR or LL (case insensitive).\n"
			"   Multiple values can be separated with commas, e.g.: 'xx,xy,yx,yy'. Two or four polarizations can be\n"
			"   joinedly cleaned (see '-joinpolarizations'), but this is not the default. I, Q, U and V\n"
			"   polarizations will be directly calculated from the visibilities, which is not appropriate for\n"
			"   telescopes with non-orthogonal feeds, such as MWA and LOFAR. The 'xy' polarization will output both\n"
			"   a real and an imaginary image, which allows calculating true Stokes polarizations for those\n"
			"   telescopes.\n"
			"-interval <start-index> <end-index>\n"
			"   Only image the given time interval. Indices specify the timesteps, end index is exclusive.\n"
			"   Default: image all time steps.\n"
			"-channelrange <start-channel> <end-channel>\n"
			"   Only image the given channel range. Indices specify channel indices, end index is exclusive.\n"
			"   Default: image all channels.\n"
			"-field <fieldid>\n"
			"   Image the given field id. Default: first field (id 0).\n"
			"-imaginarypart\n"
			"   saves the imaginary part instead of the real part; only sensible for xy/yx. Not the default.\n"
			"-datacolumn <columnname>\n"
			"   Default: CORRECTED_DATA if it exists, otherwise DATA will be used.\n"
			"-maxuvw <meters>\n"
			"-minuvw <meters>\n"
			"   Set the min/max baseline distance in meters.\n"
			"-maxw <percentage>\n"
			"   Do not grid visibilities with a w-value higher than the given percentage of the max w, to save speed.\n"
			"   Default: grid everything\n"
			"\n"
			"  ** CLEANING OPTIONS **\n"
			"-niter <niter>\n"
			"   Maximum number of clean iterations to perform. Default: 0\n"
			"-threshold <threshold>\n"
			"   Stopping clean thresholding in Jy. Default: 0.0\n"
			"-gain <gain>\n"
			"   Cleaning gain: Ratio of peak that will be subtracted in each iteration. Default: 0.1\n"
			"-mgain <gain>\n"
			"   Cleaning gain for major iterations: Ratio of peak that will be subtracted in each major\n"
			"   iteration. To use major iterations, 0.85 is a good value. Default: 1.0\n"
			"-joinpolarizations\n"
			"   Perform cleaning by searching for peaks in the sum of squares of the polarizations (either\n"
			"   I^2+Q^2+U^2+V^2, XX^2+real(XY)^2+imag(XY)^2+YY^2, or XX^2+YY^2), but subtract components from\n"
			"   individual channels. Only possible when imaging all Stokes or all linear parameters. Default: off.\n"
			"-joinchannels\n"
			"   Perform cleaning by searching for peaks in the MFS image, but subtract components from individual channels.\n"
			"   This will turn on mfsweighting by default. Default: off.\n"
			"-multiscale\n"
			"   Clean on different scales. This is a new experimental algorithm. Default: off.\n"
			"-multiscale-threshold-bias\n"
			"   Parameter to lower the threshold for larger scales. The used threshold for a scale\n"
			"   is threshold(scale)=pointsource_threshold x tbias^scale. A lower bias will clean\n"
			"   larger scales deeper. Default: 0.7\n"
			"-multiscale-scale-bias\n"
			"   Parameter to prevent cleaning small scales in the large-scale iterations. A higher\n"
			"   bias will give more focus to larger scales. Default: 0.6\n"
			"-cleanborder <percentage>\n"
			"   Set the border size in which no cleaning is performed, in percentage of the width/height of the image.\n"
			"   With an image size of 1000 and clean border of 1%, each border is 10 pixels. Default: 5 (%).\n"
			"-smallpsf\n"
			"   Resize the psf to speed up minor clean iterations. Not the default.\n"
			"-nonegative\n"
			"   Do not allow negative components during cleaning. Not the default.\n"
			"-negative\n"
			"   Default on: opposite of -nonegative.\n"
			"-stopnegative\n"
			"   Stop on negative components. Not the default.\n"
			"\n"
			"  ** RESTORATION OPTIONS **\n"
			"-beamsize <arcmin>\n"
			"   Set the FWHM beam size in arcmin for restoring the clean components. Default: longest projected\n"
			"   baseline defines restoring beam.\n";
		return -1;
	}
	
	WSClean wsclean;
	int argi = 1;
	bool mfsWeighting = false, noMFSWeighting = false, predictionMode = false;
	while(argi < argc && argv[argi][0] == '-')
	{
		const std::string param = argv[argi][1]=='-' ? (&argv[argi][2]) : (&argv[argi][1]);
		if(param == "version")
		{
			// version already printed: just exit
			return 0;
		}
		else if(param == "predict")
		{
			predictionMode = true;
		}
		else if(param == "size")
		{
			size_t
				width = atoi(argv[argi+1]),
				height = atoi(argv[argi+2]);
			if(width != height)
				throw std::runtime_error("width != height : Can't handle non-square images yet");
			wsclean.SetImageSize(width, height);
			argi += 2;
		}
		else if(param == "scale")
		{
			++argi;
			wsclean.SetPixelScale(atof(argv[argi]) * M_PI / 180.0);
		}
		else if(param == "nwlayers")
		{
			++argi;
			wsclean.SetNWlayers(atoi(argv[argi]));
		}
		else if(param == "gain")
		{
			++argi;
			wsclean.SetCleanGain(atof(argv[argi]));
		}
		else if(param == "mgain")
		{
			++argi;
			wsclean.SetCleanMGain(atof(argv[argi]));
		}
		else if(param == "niter")
		{
			++argi;
			wsclean.SetNIter(atoi(argv[argi]));
		}
		else if(param == "threshold")
		{
			++argi;
			wsclean.SetThreshold(atof(argv[argi]));
		}
		else if(param == "datacolumn")
		{
			++argi;
			wsclean.SetColumnName(argv[argi]);
		}
		else if(param == "pol")
		{
			++argi;
			wsclean.SetPolarizations(Polarization::ParseList(argv[argi]));
		}
		else if(param == "imaginarypart")
		{
			throw std::runtime_error("-imaginarypart is deprecated: imaging xy/yx will always make imaginary part.");
		}
		else if(param == "negative")
		{
			wsclean.SetAllowNegative(true);
		}
		else if(param == "nonegative")
		{
			wsclean.SetAllowNegative(false);
		}
		else if(param == "stopnegative")
		{
			wsclean.SetStopOnNegative(true);
		}
		else if(param == "makepsf")
		{
			wsclean.SetMakePSF(true);
		}
		else if(param == "cleanareas")
		{
			++argi;
			wsclean.SetCleanAreasFilename(argv[argi]);
		}
		else if(param == "name")
		{
			++argi;
			wsclean.SetPrefixName(argv[argi]);
		}
		else if(param == "gridmode")
		{
			++argi;
			std::string gridModeStr = argv[argi];
			boost::to_lower(gridModeStr);
			if(gridModeStr == "kb" || gridModeStr == "kaiserbessel" || gridModeStr == "kaiser-bessel")
				wsclean.SetGridMode(LayeredImager::KaiserBessel);
			else if(gridModeStr == "nn" || gridModeStr == "nearestneighbour")
				wsclean.SetGridMode(LayeredImager::NearestNeighbour);
			else
				throw std::runtime_error("Invalid gridding mode: should be either kb (Kaiser-Bessel) or nn (NearestNeighbour)");
		}
		else if(param == "smallinversion")
		{
			wsclean.SetSmallInversion(true);
		}
		else if(param == "nosmallinversion")
		{
			wsclean.SetSmallInversion(false);
		}
		else if(param == "smallpsf")
		{
			wsclean.SetSmallPSF(true);
		}
		else if(param == "interval")
		{
			wsclean.SetIntervalSelection(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "channelrange")
		{
			wsclean.SetChannelSelection(atoi(argv[argi+1]), atoi(argv[argi+2]));
			argi += 2;
		}
		else if(param == "channelsout")
		{
			++argi;
			wsclean.SetChannelsOut(atoi(argv[argi]));
		}
		else if(param == "joinpolarizations")
		{
			wsclean.SetJoinPolarizations(true);
		}
		else if(param == "joinchannels")
		{
			wsclean.SetJoinChannels(true);
		}
		else if(param == "mfsweighting")
		{
			mfsWeighting = true;
		}
		else if(param == "multiscale")
		{
			wsclean.SetMultiscale(true);
		}
		else if(param == "multiscale-threshold-bias")
		{
			++argi;
			wsclean.SetMultiscaleThresholdBias(atof(argv[argi]));
		}
		else if(param == "multiscale-scale-bias")
		{
			++argi;
			wsclean.SetMultiscaleScaleBias(atof(argv[argi]));
		}
		else if(param == "cleanborder")
		{
			++argi;
			wsclean.SetCleanBorderRatio(atof(argv[argi])*0.01);
		}
		else if(param == "nomfsweighting")
		{
			noMFSWeighting = true;
		}
		else if(param == "joinchannels")
		{
			wsclean.SetJoinChannels(true);
		}
		else if(param == "field")
		{
			++argi;
			wsclean.SetFieldSelection(atoi(argv[argi]));
		}
		else if(param == "weight")
		{
			++argi;
			std::string weightArg = argv[argi];
			if(weightArg == "natural")
				wsclean.SetWeightMode(WeightMode::NaturalWeighted);
			else if(weightArg == "mwa")
				wsclean.SetWeightMode(WeightMode::DistanceWeighted);
			else if(weightArg == "uniform")
				wsclean.SetWeightMode(WeightMode::UniformWeighted);
			else if(weightArg == "briggs")
			{
				++argi;
				wsclean.SetBriggsWeighting(atof(argv[argi]));
			}
			else throw std::runtime_error("Unknown weighting mode specified");
		}
		else if(param == "superweight")
		{
			++argi;
			wsclean.SetSuperWeight(atof(argv[argi]));
		}
		else if(param == "beamsize")
		{
			++argi;
			wsclean.SetBeamSize(atof(argv[argi]));
		}
		else if(param == "gkernelsize")
		{
			++argi;
			wsclean.SetAntialiasingKernelSize(atoi(argv[argi]));
		}
		else if(param == "oversampling")
		{
			++argi;
			wsclean.SetOversamplingFactor(atoi(argv[argi]));
		}
		else if(param == "reorder")
		{
			wsclean.SetForceReorder(true);
			wsclean.SetForceNoReorder(false);
		}
		else if(param == "no-reorder")
		{
			wsclean.SetForceNoReorder(true);
			wsclean.SetForceReorder(false);
		}
		else if(param == "j")
		{
			++argi;
			wsclean.SetThreadCount(atoi(argv[argi]));
		}
		else if(param == "mem")
		{
			++argi;
			wsclean.SetMemFraction(atof(argv[argi]) / 100.0);
		}
		else if(param == "absmem")
		{
			++argi;
			wsclean.SetMemAbsLimit(atof(argv[argi]));
		}
		else if(param == "minuvw")
		{
			++argi;
			wsclean.SetMinUVW(atof(argv[argi]));
		}
		else if(param == "maxuvw")
		{
			++argi;
			wsclean.SetMaxUVW(atof(argv[argi]));
		}
		else if(param == "maxw")
		{
			// This was to test the optimization suggested in Tasse et al., 2013, Appendix C.
			++argi;
			wsclean.SetWLimit(atof(argv[argi]));
		}
		else {
			throw std::runtime_error("Unknown parameter: " + param);
		}
		
		++argi;
	}
	
	if(argi == argc)
		throw std::runtime_error("No input measurement sets given.");
	
	wsclean.SetMFSWeighting((wsclean.JoinedFrequencyCleaning() && !noMFSWeighting) || mfsWeighting);
	
	for(int i=argi; i != argc; ++i)
		wsclean.AddInputMS(argv[i]);
	
	std::ostringstream commandLineStr;
	commandLineStr << "wsclean";
	for(int i=1; i!=argc; ++i)
		commandLineStr << ' ' << argv[i];
	wsclean.SetCommandLine(commandLineStr.str());
	
	if(predictionMode)
		wsclean.RunPredict();
	else
		wsclean.RunClean();
}
