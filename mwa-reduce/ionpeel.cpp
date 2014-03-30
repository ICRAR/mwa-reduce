#include <iostream>

#include "ionpeeler.h"

/**
 * Approach:
 * - Read a few timesteps (the solution interval) from the measurement set
 * - For each source in the model (multi-thread this over channels)
 *   + Predict its visibility (using Predicter)
 *   + Find the global ionospheric phase term
 *   + Subtract the source with i-term from the visibilities
 * - Write the residual flux back to the measurement set
 */

int main(int argc, char* argv[])
{
	if(argc < 4) {
		std::cout <<
			"syntax: ionpeel [options] <ms> <model> <solutionfile>\n"
			"Options are:\n"
			" -noapplybeam\n"
			"\tDon't apply the beam to the model\n"
			" -datacolumn <column>\n"
			"\tPeel from the specific column\n"
			" -weight <gridsize> <pixelscale> <mode> [robustness]\n"
			"\tApply image weights (syntax equal to wsclean)\n"
			" -groupchannels <count>\n"
			"\tGroup channels together during solving.\n"
			" -climit <flux value>\n"
			"\tRemove clusters with less than \"flux\" values.\n"
			" -distlimit <radius in deg>\n"
			"\tRemove clusters further away than given radius.\n"
			" -savemodel <out-filename>\n"
			"\tSave the model after heuristics have been applied.\n"
			" -v\n"
			"\tBe verbose.\n";
		return -1;
	}
	
	int argi = 1;
	bool applyBeam = true;
	std::string dataColumnName, outModelFilename;
	size_t solutionInterval = 1;
	WeightMode weightMode(WeightMode::NaturalWeighted);
	size_t weightGridSize = 0;
	double weightPixelScale = 0.0;
	size_t channelBlockSize = 1;
	double clusterFluxLimit = 0.0, distLimit = 0.0;
	bool verbose = false;
	
	while(argv[argi][0] == '-')
	{
		std::string param = &argv[argi][1];
		if(param == "applybeam")
		{
			applyBeam = true;
		}
		else if(param == "noapplybeam")
		{
			applyBeam = false;
		}
		else if(param == "datacolumn")
		{
			++argi;
			dataColumnName = argv[argi];
		}
		else if(param == "t")
		{
			++argi;
			solutionInterval = atoi(argv[argi]);
		}
		else if(param == "v")
		{
			verbose = true;
		}
		else if(param == "weight")
		{
			++argi;
			weightGridSize = atoi(argv[argi]);
			++argi;
			weightPixelScale = atof(argv[argi]);
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
		else if(param == "groupchannels")
		{
			++argi;
			channelBlockSize = atoi(argv[argi]);
		}
		else if(param == "climit")
		{
			++argi;
			clusterFluxLimit = atof(argv[argi]);
		}
		else if(param == "distlimit")
		{
			++argi;
			distLimit = atof(argv[argi]);
		}
		else if(param == "savemodel")
		{
			++argi;
			outModelFilename = argv[argi];
		}
		else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
		++argi;
	}
	
	if(clusterFluxLimit != 0.0 && outModelFilename.empty())
	{
		throw std::runtime_error("You have not specified an output model filename (with -savemodel) but have specified a cluster flux limit (-climit): you will not be able to apply the solutions without the changed model!");
	}

	IonPeeler peeler;
	peeler.SetSolutionInterval(solutionInterval);
	peeler.SetApplyBeam(applyBeam);
	peeler.SetDataColumnName(dataColumnName);
	peeler.SetWeighting(weightMode, weightGridSize, weightPixelScale);
	peeler.SetChannelBlockSize(channelBlockSize);
	peeler.SetVerbose(verbose);
	if(clusterFluxLimit != 0.0)
		peeler.SetClusterFluxLimit(clusterFluxLimit);
	if(distLimit != 0.0)
		peeler.SetDistanceLimit(distLimit);
	peeler.Peel(argv[argi], argv[argi+1], argv[argi+2]);
	
	if(!outModelFilename.empty())
		peeler.GetUsedModel().Save(outModelFilename.c_str());
}
