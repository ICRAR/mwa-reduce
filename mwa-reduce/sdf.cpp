#include <iostream>

#include "banddata.h"
#include "model.h"

#include <ms/MeasurementSets/MeasurementSet.h>

int main(int argc, char *argv[])
{
	int argi = 1;
	bool outputPlot = false;
	if(strcmp(argv[argi], "-p") == 0)
	{
		outputPlot = true;
		++argi;
	}
	size_t newChannelCount = atoi(argv[argi]);
	Model model(argv[argi+1]);
	casa::MeasurementSet ms(argv[argi+2]);
	BandData band(ms.spectralWindow());
	const long double newBandSize = (band.BandEnd() - band.BandStart()) / newChannelCount;
	
	for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
	{
		SourceSDFWithSamples<long double> newSdf;
		const SourceSDF<long double> &sdf = sourcePtr->Brightness();
		for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
		{
			long double startFreq = band.BandStart() + newBandSize*newChIndex;
			long double endFreq = band.BandStart() + newBandSize*(newChIndex+1.0);
			long double flux = sdf.IntegratedFlux(startFreq, endFreq);
			newSdf.AddSample(flux, (startFreq+endFreq)*0.5);
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
