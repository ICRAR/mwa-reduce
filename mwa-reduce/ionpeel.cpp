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
	if(argc < 3) {
		std::cout << "syntax: ionpeel [options] <ms> <model>\n";
		return -1;
	}
	
	int argi = 1;
	bool applyBeam = true;
	std::string dataColumnName = "DATA";
	size_t solutionInterval = 1;
	WeightMode weightMode(WeightMode::NaturalWeighted);
	size_t weightGridSize = 0;
	double weightPixelScale = 0.0;
	
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
		else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
		++argi;
	}

	IonPeeler peeler;
	peeler.SetSolutionInterval(solutionInterval);
	peeler.SetApplyBeam(applyBeam);
	peeler.SetDataColumnName(dataColumnName);
	peeler.SetWeighting(weightMode, weightGridSize, weightPixelScale);
	peeler.Peel(argv[argi], argv[argi+1]);
}
