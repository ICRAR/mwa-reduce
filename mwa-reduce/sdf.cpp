#include "banddata.h"
#include "model.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		std::cout << "sdf -- Interpolation, extrapolation, plotting and scaling of the \n"
		"spectral density function. Usage:\n"
		"\tsdf [-p] [-m <output model>] [-o] [-s <scale>] [-set0/1/2/3 <flux>] [-unpolarized] [-pl] [-t <threshold>] [-r <new-nr-channels>] <model>\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, powerlaw = false, optimize = false, applyThreshold = false, resample = false, unpolarized = false;
	bool setPolarization[4] = {false, false, false, false};
	long double setPolFlux[4] = {0.0, 0.0, 0.0, 0.0};
	long double scale = 1.0, threshold = 0.0;
	size_t newChannelCount = 0;
	std::string outputModel;
	while(argv[argi][0]=='-')
	{
		if(strcmp(argv[argi], "-p") == 0)
		{
			outputPlot = true;
		} else if(strcmp(argv[argi], "-m") == 0)
		{
			++argi;
			outputModel = argv[argi];
		} else if(strcmp(argv[argi], "-s") == 0)
		{
			++argi;
			scale = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-set0") == 0)
		{
			++argi;
			setPolarization[0] = true;
			setPolFlux[0] = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-set1") == 0)
		{
			++argi;
			setPolarization[1] = true;
			setPolFlux[1] = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-set2") == 0)
		{
			++argi;
			setPolarization[2] = true;
			setPolFlux[2] = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-set3") == 0)
		{
			++argi;
			setPolarization[3] = true;
			setPolFlux[3] = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-pl") == 0)
		{
			powerlaw = true;
		} else if(strcmp(argv[argi], "-t") == 0)
		{
			++argi;
			threshold = atof(argv[argi]);
			applyThreshold = true;
		} else if(strcmp(argv[argi], "-r") == 0)
		{
			++argi;
			newChannelCount = atoi(argv[argi]);
			resample = true;
		} else if(strcmp(argv[argi], "-unpolarized") == 0)
		{
			unpolarized = true;
		} else if(strcmp(argv[argi], "-o") == 0)
		{
			optimize = true;
		} else {
			throw std::runtime_error("Unknown option given");
		}
		++argi;
	}
	Model model(argv[argi]);
	
	if(optimize)
		model.Optimize();
	
	if(applyThreshold)
	{
		Model temp;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const SpectralEnergyDistribution &sed = sourcePtr->SED();
			if(sed.IntegratedFlux(sed.LowestFrequency(), sed.HighestFrequency(), 0) >= threshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
	{
		const SpectralEnergyDistribution &sed = sourcePtr->SED();
		const long double startFreq = sed.LowestFrequency();
		const long double endFreq = sed.HighestFrequency();
		if(powerlaw)
		{
			long double e, f;
			SpectralEnergyDistribution newSED;
			Measurement m1, m2;
			m1.SetFrequencyHz(startFreq);
			m2.SetFrequencyHz(endFreq);
			
			for(size_t p=0; p!=4; ++p)
			{
				sed.FitPowerlaw(f, e, p);
				m1.SetFluxDensity(p, f * std::pow(startFreq, e));
				m2.SetFluxDensity(p, f * std::pow(endFreq, e));
			}
			newSED.AddMeasurement(m1);
			newSED.AddMeasurement(m2);
			sourcePtr->SetSED(newSED);
		}
		else if(resample)
		{
			SpectralEnergyDistribution newSED;
			const long double newBandSize = (endFreq - startFreq) / newChannelCount;
			for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
			{
				long double chStartFreq = startFreq + newBandSize*newChIndex;
				long double chEndFreq = endFreq + newBandSize*(newChIndex+1.0);
				long double flux;
				if(newChannelCount < sed.MeasurementCount())
					flux = sed.IntegratedFlux(chStartFreq, chEndFreq, 0);
				else
					flux = sed.FluxAtFrequency((chStartFreq+chEndFreq)*0.5, 0);
				
				newSED.AddMeasurement(flux*scale, (chStartFreq+chEndFreq)*0.5);
			}
			
			sourcePtr->SetSED(newSED);
		}
	}
	
	for(size_t p=0; p!=4; ++p)
	{
		if(setPolarization[p])
		{
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				SpectralEnergyDistribution &sed = sourcePtr->SED();
				for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
				{
					m->SetFluxDensity(p, setPolFlux[p]);
				}
			}
		}
	}
	
	if(unpolarized)
	{
		model.SetUnpolarized();
	}
	
	if(outputPlot)
	{
		std::ofstream plotStream("spectrum.plt");
		plotStream <<
			"set terminal postscript enhanced color\n"
			"#set logscale xy\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"spectrum.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n"
			"plot \\\n";

		std::ofstream plotIStream("spectrum-I.plt");
		plotIStream <<
			"set terminal postscript enhanced color\n"
			"#set logscale xy\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"spectrum-I.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n"
			"plot \\\n";

		for(size_t sourceIndex = 0; sourceIndex!=model.SourceCount(); ++sourceIndex)
		{
			std::ostringstream dataStreamName;
			dataStreamName << "spectrum" << sourceIndex << ".txt";
			std::ofstream dataStream(dataStreamName.str().c_str());
			plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 2.0 title \"\",\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:3 with lines lw 2.0 title \"\",\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:4 with lines lw 2.0 title \"\",\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:5 with lines lw 2.0 title \"\"";
			
			plotIStream << "\"" << dataStreamName.str() << "\" using 1:((column(2)+column(5))*0.5) with lines lw 2.0 title \"\",\\\n";
			
			const SpectralEnergyDistribution &sed = model.Source(sourceIndex).SED();
			long double e1, e2, f1, f2;
			sed.FitPowerlaw(f1, e1, 0);
			sed.FitPowerlaw(f2, e2, 3);
			plotIStream << (f1/2.0) << " * (x*1000000)**" << e1 << " + " << (f2/2.0) << " * (x*1000000)**" << e2 << " with lines lw 1.0 title \"\"";
			
			if(sourceIndex != model.SourceCount()-1)
			{
				plotStream << ",\\";
				plotIStream << ",\\";
			}
			plotStream << "\n";
			plotIStream << "\n";
			
			std::vector<Measurement> measurements;
			sed.GetMeasurements(measurements);
			
			for(std::vector<Measurement>::const_iterator i=measurements.begin(); i!=measurements.end(); ++i)
			{
				dataStream
					<< i->FrequencyHz()/1000000.0 << '\t'
					<< i->FluxDensity(0) << '\t'
					<< i->FluxDensity(1) << '\t'
					<< i->FluxDensity(2) << '\t'
					<< i->FluxDensity(3) << '\n';
			}
		}
	}
	
	if(!outputModel.empty())
		model.Save(outputModel.c_str());
	
	return 0;
}
