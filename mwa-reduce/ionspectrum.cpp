#include "ionspectrummaker.h"

int main(int argc, char* argv[])
{
	if(argc < 4) {
		std::cout << "syntax: ionspectrum [options] <ms> <ionsolutions> <model> <output>\n"
			"Options:\n"
			" -weight <gridsize> <pixelscale> <mode> [robustness]\n";
		return -1;
	}
	
	WeightMode weightMode(WeightMode::NaturalWeighted);
	size_t weightGridSize = 0;
	double weightPixelScale = 0.0;
	int argi = 1;
	while(argv[argi][0] == '-')
	{
		std::string param = &argv[argi][1];
		if(param == "weight")
		{
			++argi;
			weightGridSize = atoi(argv[argi]);
			++argi;
			weightPixelScale = atof(argv[argi]) * M_PI / 180.0;
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
		else throw std::runtime_error("Unknown option specified");
	}	
	const char *msFilename(argv[argi]);
	const char *ionFilename(argv[argi+1]);
	const char *modelFilename(argv[argi+2]);
	const char *outputFilename(argv[argi+3]);
	IonSpectrumMaker isMaker;
	isMaker.SetWeighting(weightMode, weightGridSize, weightPixelScale);
	isMaker.InitializeForVisibilityAcc(msFilename, ionFilename, modelFilename);
	isMaker.Accumulate();
	isMaker.Save(outputFilename);
	
	return 0;
}
