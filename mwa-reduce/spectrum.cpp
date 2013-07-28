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
		std::cout << "Usage: spectrum <model> <ms>\n"
			"Calculates the spectrum directly from the ms, for each source in the model.\n";
	} else {
		size_t argi = 1;
		Model model(argv[argi]);
		Model measuredModel(model);
		
		SpectrumMaker spectrumMaker;
		spectrumMaker.AddMeasurementSet(argv[argi+1]);
		
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
			spectrumMaker.AddSource(*s);
		
		spectrumMaker.Measure();
		
		bool outputModel = true;
		if(outputModel)
		{
			size_t sourceIndex = 0;
			for(Model::iterator source=model.begin();source!=model.end();++source)
			{
				SpectralEnergyDistribution sed;
				size_t itemIndex = sourceIndex;
				for(size_t ch=0; ch!=channelCount;++ch)
				{
					sed.AddMeasurement(sourceFlux[itemIndex] / sourceMeasCount[itemIndex], bandData.ChannelFrequency(ch));
					itemIndex += model.SourceCount();
				}
				source->SetSED(sed);
				std::cout << source->ToString() << '\n';
				++sourceIndex;
			}
		} else {
			float sums[model.SourceCount()];
			size_t counts[model.SourceCount()];
			for(size_t i=0;i!=model.SourceCount();++i)
			{
				sums[i] = 0;
				counts[i] = 0;
			}
			
			for(size_t ch=0; ch!=channelCount;++ch)
			{
				std::cout << ch << '\t' << bandData.ChannelFrequency(ch);
				size_t sourceIndex = ch * model.SourceCount();
				for(size_t s=0; s!=model.SourceCount();++s)
				{
					std::cout << '\t' << (sourceFlux[sourceIndex] / sourceMeasCount[sourceIndex]);
					sums[s] += sourceFlux[sourceIndex];
					counts[s] += sourceMeasCount[sourceIndex];
					++sourceIndex;
				}
				std::cout << '\n';
			}
			std::cout << "avg\tavg\t";
			for(size_t i=0;i!=model.SourceCount();++i)
			{
				std::cout << '\t' << (sums[i] / counts[i]);
			}
			std::cout << '\n';
		}
	}
}
