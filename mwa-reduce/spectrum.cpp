#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "banddata.h"
#include "sourcesdf.h"
#include "model.h"
#include "predicter.h"
#include "spectrummaker.h"

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: spectrum [-applybeam] [-s <model to subtract>] [-g <solutions>] [-saveintermediate <file>] <model for positions> <output-model> <ms> [<ms2>...]\n"
			"Calculates the spectrum directly from the ms, for each source in the model.\n";
	} else {
		size_t argi = 1;
		const char *subtractionModelFile = 0, *solutionsFile = 0, *intermediateFile = 0;
		bool applyBeam = false;
		while(argv[argi][0] == '-')
		{
			const std::string param = &argv[argi][1];
			if(param == "s")
			{
				++argi;
				subtractionModelFile = argv[argi];
			}
			else if(param == "g")
			{
				++argi;
				solutionsFile = argv[argi];
			}
			else if(param == "applybeam")
			{
				applyBeam = true;
			}
			else if(param == "saveintermediate")
			{
				++argi;
				intermediateFile = argv[argi];
			}
			else throw std::runtime_error("Invalid parameter");
			++argi;
		}
		
		Model model(argv[argi]);
		const char *modelFilename = argv[argi+1];
		
		SpectrumMaker spectrumMaker;
		for(int i=argi+2; i!=argc; ++i)
		{
			if(solutionsFile != 0)
				spectrumMaker.AddMeasurementSet(argv[i], solutionsFile);
			else
				spectrumMaker.AddMeasurementSet(argv[i]);
		}
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
			spectrumMaker.AddSource(*s);
		spectrumMaker.SetApplyBeam(applyBeam);
		
		if(subtractionModelFile != 0) spectrumMaker.SetSubtractedModel(Model(subtractionModelFile));
		
		spectrumMaker.Measure();
		
		if(intermediateFile != 0)
			spectrumMaker.SaveIntermediate(intermediateFile);
		
		spectrumMaker.Finish();
		
		Model outputModel;
		spectrumMaker.ToModel(outputModel);
		outputModel.Save(modelFilename);
	}
}
