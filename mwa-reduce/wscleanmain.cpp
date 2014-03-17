#include "wsclean.h"

#include <boost/algorithm/string.hpp>

#include <iostream>

int main(int argc, char *argv[])
{
	std::cout << "\n"
		" ** This software package is released under the GPL version 3. **\n"
	  " ** Author: Andre Offringa (offringa@gmail.com).               **\n";
	
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
			"\t-smallinversion\n"
			"\t   Experimental mode to speed up inversion.\n"
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
			"\t   Default: grid everything\n"
			"\t-mem <percentage>\n"
			"\t   Limit memory usage to the given fraction of the total system memory. This is an approximate value.\n"
			"\t   Default: 100.\n";
		return -1;
	}
	
	WSClean wsclean;
	int argi = 1;
	while(argi < argc && argv[argi][0] == '-')
	{
		const std::string param = &argv[argi][1];
		if(param == "size")
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
		else if(param == "addmodel")
		{
			++argi;
			wsclean.SetAddModelFilename(argv[argi]);
		}
		else if(param == "addmodelapp")
		{
			++argi;
			wsclean.SetAddModelFilename(argv[argi]);
			wsclean.SetAddAppModel(true);
		}
		else if(param == "savemodel")
		{
			++argi;
			wsclean.SetSaveModelFilename(argv[argi]);
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
		else if(param == "mem")
		{
			++argi;
			wsclean.SetMemFraction(atof(argv[argi]) / 100.0);
		}
		else if(param == "wlimit")
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
	
	for(int i=argi; i != argc; ++i)
		wsclean.AddInputMS(argv[i]);
	
	std::ostringstream commandLineStr;
	commandLineStr << "wsclean";
	for(int i=1; i!=argc; ++i)
		commandLineStr << ' ' << argv[i];
	wsclean.SetCommandLine(commandLineStr.str());
			
	wsclean.Run();
}
