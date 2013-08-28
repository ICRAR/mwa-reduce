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
		std::cout << "Usage: spectrum [-s <model to subtract>] [-g <solutions>] [-applybeam] <model for positions> <output-model> <ms> [<ms2>...]\n"
			"Calculates the spectrum directly from the ms, for each source in the model.\n";
	} else {
		size_t argi = 1;
		const char *subtractionModelFile = 0, *solutionsFile = 0;
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
		
		Model outputModel;
		for(size_t sourceIndex = 0; sourceIndex!=model.SourceCount(); ++sourceIndex)
		{
			SpectralEnergyDistribution sed;
			std::map<double, double> spectrum[4];
			for(size_t p=0; p!=4; ++p)
				spectrumMaker.FluxPerFrequency(spectrum[p], sourceIndex, p);
			
			std::map<double, double>::const_iterator
				chIter1 = spectrum[1].begin(), chIter2 = spectrum[2].begin(), chIter3 = spectrum[3].begin();
			for(std::map<double, double>::const_iterator chIter0=spectrum[0].begin(); chIter0!=spectrum[0].end(); ++chIter0)
			{
				Measurement m;
				m.SetFluxDensity(0, chIter0->second);
				m.SetFluxDensity(1, chIter1->second);
				m.SetFluxDensity(2, chIter2->second);
				m.SetFluxDensity(3, chIter3->second);
				m.SetFrequencyHz(chIter0->first);
				sed.AddMeasurement(m);
				
				++chIter1; ++chIter2; ++chIter3;
			}
			ModelSource source = model.Source(sourceIndex);
			source.SetSED(sed);
			outputModel.AddSource(source);
		}
		outputModel.Save(modelFilename);
	}
}
