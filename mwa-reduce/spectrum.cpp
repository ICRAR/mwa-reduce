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
		std::cout << "Usage: spectrum <model> <output-model> <ms> [<ms2>...]\n"
			"Calculates the spectrum directly from the ms, for each source in the model.\n";
	} else {
		size_t argi = 1;
		Model model(argv[argi]);
		const char *modelFilename = argv[argi+1];
		
		SpectrumMaker spectrumMaker;
		for(int i=argi+2; i!=argc; ++i)
		{
			spectrumMaker.AddMeasurementSet(argv[i]);
		}
		
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
			spectrumMaker.AddSource(*s);
		
		spectrumMaker.Measure();
		
		Model outputModel;
		for(size_t sourceIndex = 0; sourceIndex!=model.SourceCount(); ++sourceIndex)
		{
			SpectralEnergyDistribution sed;
			std::map<double, long double> spectrum[4];
			for(size_t p=0; p!=4; ++p)
				spectrumMaker.FluxPerFrequency(spectrum[p], sourceIndex, p);
			
			std::map<double, long double>::const_iterator
				chIter1 = spectrum[1].begin(), chIter2 = spectrum[2].begin(), chIter3 = spectrum[3].begin();
			for(std::map<double, long double>::const_iterator chIter0=spectrum[0].begin(); chIter0!=spectrum[0].end(); ++chIter0)
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
		
		std::ofstream plotStream("spectrum.plt");
		plotStream <<
			"set terminal postscript enhanced color\n"
			"#set logscale y\n"
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
			"#set logscale y\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"spectrum.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n"
			"plot \\\n";

		for(size_t sourceIndex = 0; sourceIndex!=model.SourceCount(); ++sourceIndex)
		{
			std::ostringstream dataStreamName;
			dataStreamName << "spectrum" << sourceIndex << ".txt";
			std::ofstream dataStream(dataStreamName.str().c_str());
			plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 2.0,\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:3 with lines lw 2.0,\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:4 with lines lw 2.0,\\\n";
			plotStream << "\"" << dataStreamName.str() << "\" using 1:5 with lines lw 2.0";
			plotIStream << "\"" << dataStreamName.str() << "\" using 1:((column(2)+column(4))*0.5) with lines lw 2.0";
			if(sourceIndex != model.SourceCount()-1)
			{
				plotStream << ",\\";
				plotIStream << ",\\";
			}
			plotStream << "\n";
			plotIStream << "\n";
			std::map<double, long double> spectrum[4];
			for(size_t p=0; p!=4; ++p)
				spectrumMaker.FluxPerFrequency(spectrum[p], sourceIndex, p);
			std::map<double, long double>::const_iterator
				chIter1 = spectrum[1].begin(), chIter2 = spectrum[2].begin(), chIter3 = spectrum[3].begin();
				
			for(std::map<double, long double>::const_iterator chIter0=spectrum[0].begin(); chIter0!=spectrum[0].end(); ++chIter0)
			{
				dataStream
					<< chIter0->first/1000000.0 << '\t'
					<< chIter0->second << '\t'
					<< chIter1->second << '\t'
					<< chIter2->second << '\t'
					<< chIter3->second << '\n';
				++chIter1; ++chIter2; ++chIter3;
			}
		}
	}
}
