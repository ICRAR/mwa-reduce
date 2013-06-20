#include <iostream>

#include "banddata.h"
#include "model.h"

#include <ms/MeasurementSets/MeasurementSet.h>

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		std::cout << "sdf -- Interpolation, extrapolation, plotting and scaling of the \n"
		"spectral density function. Usage:\n"
		"\tsdf [-p] [-o] [-s <scale>] [-t <threshold>] <new-nr-channels> <model> <ms>\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, outputSIs = false, optimize = false, applyThreshold = false;
	long double scale = 1.0, threshold = 0.0;
	while(argv[argi][0]=='-')
	{
		if(strcmp(argv[argi], "-p") == 0)
		{
			outputPlot = true;
		} else if(strcmp(argv[argi], "-s") == 0)
		{
			++argi;
			scale = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-si") == 0)
		{
			outputSIs = true;
		} else if(strcmp(argv[argi], "-t") == 0)
		{
			++argi;
			threshold = atof(argv[argi]);
			applyThreshold = true;
		} else if(strcmp(argv[argi], "-o") == 0)
		{
			optimize = true;
		} else {
			throw std::runtime_error("Unknown option given");
		}
		++argi;
	}
	size_t newChannelCount = atoi(argv[argi]);
	Model model(argv[argi+1]);
	casa::MeasurementSet ms(argv[argi+2]);
	BandData band(ms.spectralWindow());
	
	if(optimize)
		model.Optimize();
	if(applyThreshold)
	{
		Model temp;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			if(sourcePtr->SED().IntegratedFlux(band.BandStart(), band.BandEnd(), 0) >= threshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	const long double newBandSize = (band.BandEnd() - band.BandStart()) / newChannelCount;
	for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
	{
		const SpectralEnergyDistribution &sed = sourcePtr->SED();
		if(outputSIs)
		{
			long double startFreq = band.BandStart();
			long double endFreq = band.BandEnd();
			SpectralEnergyDistribution newSED(
				sed.FluxAtFrequency(startFreq, 0), startFreq,
				sed.FluxAtFrequency(endFreq, 0), endFreq);
			if(!outputPlot) {
				sourcePtr->SetSED(newSED);
				std::cout << sourcePtr->ToString() << '\n';
			}
		}
		else {
			SpectralEnergyDistribution newSED;
			for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
			{
				long double startFreq = band.BandStart() + newBandSize*newChIndex;
				long double endFreq = band.BandStart() + newBandSize*(newChIndex+1.0);
				long double flux;
				if(newChannelCount < band.ChannelCount())
					flux = sed.IntegratedFlux(startFreq, endFreq, 0);
				else
					flux = sed.FluxAtFrequency((startFreq+endFreq)*0.5, 0);
				
				newSED.AddMeasurement(flux*scale, (startFreq+endFreq)*0.5);
			}
			if(!outputPlot) {
				sourcePtr->SetSED(newSED);
				std::cout << sourcePtr->ToString() << '\n';
			}
		}
	}
	
	if(outputPlot)
	{
		for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
		{
			long double startFreq = band.BandStart() + newBandSize*newChIndex;
			long double endFreq = band.BandStart() + newBandSize*(newChIndex+1.0);
			std::cout << (endFreq+startFreq)*0.5;
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				const SpectralEnergyDistribution &sed = sourcePtr->SED();
				long double flux = sed.IntegratedFlux(startFreq, endFreq, 0);
				std::cout << '\t' << flux;
			}
			std::cout << '\n';
		}
	}
	
	return 0;
}
