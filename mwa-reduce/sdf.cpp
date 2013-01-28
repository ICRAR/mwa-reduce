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
		"\tsdf [-p] [-o] [-s <scale>] <new-nr-channels> <model> <ms>\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, optimize = false;
	long double scale = 1.0;
	while(argv[argi][0]=='-')
	{
		if(strcmp(argv[argi], "-p") == 0)
		{
			outputPlot = true;
		} else if(strcmp(argv[argi], "-s") == 0)
		{
			++argi;
			scale = atof(argv[argi]);
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
	const long double newBandSize = (band.BandEnd() - band.BandStart()) / newChannelCount;
	
	if(optimize)
		model.Optimize();
	
	for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
	{
		SourceSDFWithSamples<long double> newSdf;
		const SourceSDF<long double> &sdf = sourcePtr->Brightness();
		for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
		{
			long double startFreq = band.BandStart() + newBandSize*newChIndex;
			long double endFreq = band.BandStart() + newBandSize*(newChIndex+1.0);
			long double flux;
			if(newChannelCount < band.ChannelCount())
				flux = sdf.IntegratedFlux(startFreq, endFreq);
			else
				flux = sdf.FluxAtFrequency((startFreq+endFreq)*0.5);
			
			newSdf.AddSample(flux*scale, (startFreq+endFreq)*0.5);
		}
		if(!outputPlot) {
			sourcePtr->SetBrightness(newSdf);
			std::cout << sourcePtr->ToStringLine() << '\n';
		}
	}
	
	if(outputPlot)
	{
		for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
		{
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				const SourceSDF<long double> &sdf = sourcePtr->Brightness();
				long double startFreq = band.BandStart() + newBandSize*newChIndex;
				long double endFreq = band.BandStart() + newBandSize*(newChIndex+1.0);
				long double flux = sdf.IntegratedFlux(startFreq, endFreq);
				std::cout << (endFreq+startFreq)*0.5 << '\t' << flux << '\n';
			}
		}
	}
	
	return 0;
}
