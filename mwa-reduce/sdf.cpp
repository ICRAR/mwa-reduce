#include "banddata.h"
#include "model.h"
#include "imagecoordinates.h"
#include "areaset.h"
#include "parser/areaparser.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <fstream>
#include <set>

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		std::cout << "sdf -- Interpolation, extrapolation, plotting and scaling of the \n"
		"spectral density function. Usage:\n"
		"\tsdf [-p] [-m <output model>] [-o] [-s <scale>] [-sp <peakflux A> <freq A> <peakflux B> <freq B>] [-set0/1/2/3 <flux>] [-unpolarized] [-pl] [-t <threshold>] [-r <new-nr-channels>] [-delnoisysources <fluxlimit>] [-near <ra> <dec> <dist>] [-combine-diff-meas] <model> [<more models>..] [-collect <name>] [-sort] [-without/only <areafile>]\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, powerlaw = false, optimize = false, applyThreshold = false, resample = false, unpolarized = false, delNoisySources = false;
	bool setPolarization[4] = {false, false, false, false};
	long double setPolFlux[4] = {0.0, 0.0, 0.0, 0.0};
	long double scale = 1.0, threshold = 0.0, delNoisySourceLimit = 0.0;
	long double scalePeakA = 1.0, scaleFreqA = 0.0, scalePeakB = 1.0, scaleFreqB = 0.0;
	size_t newChannelCount = 0;
	std::string outputModel, collectName;
	bool nearFilter = false, scalePeak = false, scaleSource = false, doCollect = false, doSort = false;
	long double nearFilterRA = 0.0, nearFilterDec = 0.0, nearFilterDist = 0.0;
	enum { AddFluxes, AverageFluxes, DifferentFrequencies } combineStrategy = AddFluxes;
	bool excludeArea = false;
	std::unique_ptr<AreaSet> areaSet;
	while(argv[argi][0]=='-')
	{
		if(strcmp(argv[argi], "-p") == 0)
		{
			outputPlot = true;
		} else if(strcmp(argv[argi], "-near") == 0)
		{
			++argi;
			nearFilterRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			nearFilterDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			nearFilterDist = atof(argv[argi]) * (M_PI/180.0);
			nearFilter = true;
		} else if(strcmp(argv[argi], "-collect") == 0)
		{
			doCollect = true;
			++argi;
			collectName = argv[argi];
		} else if(strcmp(argv[argi], "-m") == 0)
		{
			++argi;
			outputModel = argv[argi];
		} else if(strcmp(argv[argi], "-delnoisysources") == 0)
		{
			++argi;
			delNoisySources = true;
			delNoisySourceLimit = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-s") == 0)
		{
			++argi;
			scale = atof(argv[argi]);
		} else if(strcmp(argv[argi], "-sp") == 0)
		{
			scalePeak = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(strcmp(argv[argi], "-sc") == 0)
		{
			scaleSource = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
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
		} else if(strcmp(argv[argi], "-combine-diff-meas") == 0)
		{
			combineStrategy = DifferentFrequencies;
		} else if(strcmp(argv[argi], "-combine-avg-meas") == 0)
		{
			combineStrategy = AverageFluxes;
		} else if(strcmp(argv[argi], "-sort") == 0)
		{
			doSort = true;
		} else if(strcmp(argv[argi], "-without") == 0 || strcmp(argv[argi], "-only") == 0)
		{
			excludeArea = strcmp(argv[argi], "-without") == 0;
			++argi;
			areaSet.reset(new AreaSet());
			AreaParser areaParser;
			areaParser.Parse(*areaSet, argv[argi]);
		} else {
			throw std::runtime_error(std::string("Unknown option given: ") + argv[argi]);
		}
		++argi;
	}
	Model model;
	switch(combineStrategy)
	{
		case AddFluxes:
		case AverageFluxes:
			for(int modelIndex=argi; modelIndex!=argc; ++modelIndex)
			{
				model += Model(argv[modelIndex]);
			}
			break;
		case DifferentFrequencies:
			for(int modelIndex=argi; modelIndex!=argc; ++modelIndex)
			{
				model.CombineMeasurements(Model(argv[modelIndex]));
			}
			break;
	}
	if(combineStrategy == AverageFluxes)
	{
		double fact = 1.0 / (argc - argi);
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				SpectralEnergyDistribution &sed = compPtr->SED();
				for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
				{
					for(size_t p=0; p!=4; ++p)
						m->second.SetFluxDensity(p, m->second.FluxDensity(p) * fact);
				}
			}
		}
	}
	
	if(optimize)
		model.Optimize();
	if(doSort)
		model.Sort();
	
	if(areaSet != 0)
	{
		std::vector<size_t> toBeRemoved;
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			const ModelSource& source = model.Source(i);
			std::vector<const SkyArea*> areas;
			areaSet->FindAreas(areas, source.Peak().PosRA(), source.Peak().PosDec());
			if((areas.empty() && !excludeArea) || (!areas.empty() && excludeArea))
				toBeRemoved.push_back(i);
		}
		for(std::vector<size_t>::const_reverse_iterator i=toBeRemoved.rbegin(); i!=toBeRemoved.rend(); ++i)
			model.RemoveSource(*i);
	}
	
	if(applyThreshold)
	{
		Model temp;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const SpectralEnergyDistribution& sed = sourcePtr->Peak().SED();
			double flux = sourcePtr->TotalFlux(sed.LowestFrequency(), sed.HighestFrequency(), 0);
			if(!std::isfinite(flux))
			{
				long double f, e;
				sed.FitPowerlaw(f, e, 0);
				double lowFlux = f * std::pow(sed.LowestFrequency(), e);
				double hiFlux = f * std::pow(sed.HighestFrequency(), e);
				flux = (lowFlux + hiFlux) * 0.5;
			}

			if(flux >= threshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	if(delNoisySources)
	{
		std::set<size_t> sourcesToDelete;
		for(size_t p=0; p!=4; ++p)
		{
			for(size_t i = 0; i!=model.SourceCount(); ++i)
			{
				SpectralEnergyDistribution &sed = model.Source(i).Peak().SED();
				for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
				{
					long double flux = m->second.FluxDensity(p);
					if(flux > delNoisySourceLimit || flux < -delNoisySourceLimit)
						sourcesToDelete.insert(i);
				}
			}
		}
		for(std::set<size_t>::reverse_iterator i=sourcesToDelete.rbegin(); i!=sourcesToDelete.rend(); ++i)
			model.RemoveSource(*i);
	}
	
	if(nearFilter)
	{
		for(size_t i = model.SourceCount(); i>0; --i)
		{
			ModelSource& source = model.Source(i-1);
			double dist = ImageCoordinates::AngularDistance(source.Peak().PosRA(), source.Peak().PosDec(), nearFilterRA, nearFilterDec);
			if(dist > nearFilterDist)
				model.RemoveSource(i-1);
		}
	}
	
	if(doCollect)
	{
		ModelSource newSource;
		newSource.SetName(collectName);
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr != sourcePtr->end(); ++compPtr)
			{
				newSource.AddComponent(*compPtr);
			}
		}
		model = Model();
		model.AddSource(newSource);
	}
	
	for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
	{
		const SpectralEnergyDistribution &sed = sourcePtr->Peak().SED();
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
			sourcePtr->Peak().SetSED(newSED);
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
				
				newSED.AddMeasurement(flux, (chStartFreq+chEndFreq)*0.5);
			}
			
			sourcePtr->Peak().SetSED(newSED);
		}
	}
	
	for(size_t p=0; p!=4; ++p)
	{
		if(setPolarization[p])
		{
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
				{
					SpectralEnergyDistribution &sed = compPtr->SED();
					for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
					{
						m->second.SetFluxDensity(p, setPolFlux[p]);
					}
				}
			}
		}
		if(scale != 1.0)
		{
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
				{
					SpectralEnergyDistribution &sed = compPtr->SED();
					for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
					{
						m->second.SetFluxDensity(p, m->second.FluxDensity(p) * scale);
					}
				}
			}
		}
	}
	
	if(scalePeak || scaleSource)
	{
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			long double factorA[4], factorB[4];
			for(size_t p=0; p!=4; ++p)
			{
				if(scalePeak)
				{
					long double
						oldFluxA = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqA, p),
						oldFluxB = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqB, p);
					factorA[p] = oldFluxA==0.0 ? 0.0 : scalePeakA / oldFluxA;
					factorB[p] = oldFluxB==0.0 ? 0.0 : scalePeakB / oldFluxB;
				} else {
					long double
						oldFluxA = sourcePtr->TotalFlux(scaleFreqA, p),
						oldFluxB = sourcePtr->TotalFlux(scaleFreqB, p);
					factorA[p] = oldFluxA==0.0 ? 0.0 : scalePeakA / oldFluxA;
					factorB[p] = oldFluxB==0.0 ? 0.0 : scalePeakB / oldFluxB;
				}
			}
			std::cout << "Scale factor for " << sourcePtr->Name() << ": " << factorA[0];
			for(size_t p=1; p!=4; ++p) std::cout << ',' << factorA[p];
			std::cout << " - " << factorB[0];
			for(size_t p=1; p!=4; ++p) std::cout << ',' << factorB[p];
			std::cout << '\n';
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				Measurement measA, measB;
				measA.SetFrequencyHz(scaleFreqA);
				measB.SetFrequencyHz(scaleFreqB);
				for(size_t p=0; p!=4; ++p)
				{
					long double oldFluxA = compPtr->SED().FluxAtFrequency(scaleFreqA, p);
					long double oldFluxB = compPtr->SED().FluxAtFrequency(scaleFreqB, p);
					measA.SetFluxDensity(p, oldFluxA*factorA[p]);
					measB.SetFluxDensity(p, oldFluxB*factorB[p]);
				}
				SpectralEnergyDistribution sed;
				sed.AddMeasurement(measA);
				sed.AddMeasurement(measB);
				compPtr->SetSED(sed);
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
			"set logscale y\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"spectrum-I.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n"
			"plot \\\n";

		size_t sourceIndex = 0;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				std::ostringstream dataStreamName;
				dataStreamName << "spectrum" << sourceIndex << ".txt";
				std::ofstream dataStream(dataStreamName.str().c_str());
				plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 2.0 title \"\",\\\n";
				plotStream << "\"" << dataStreamName.str() << "\" using 1:3 with lines lw 2.0 title \"\",\\\n";
				plotStream << "\"" << dataStreamName.str() << "\" using 1:4 with lines lw 2.0 title \"\",\\\n";
				plotStream << "\"" << dataStreamName.str() << "\" using 1:5 with lines lw 2.0 title \"\"";
				
				plotIStream << "\"" << dataStreamName.str() << "\" using 1:((column(2)+column(5))*0.5) with lines lw 2.0 title \"\",\\\n";
				
				const SpectralEnergyDistribution &sed = compPtr->SED();
				long double e1, e2, f1, f2;
				sed.FitPowerlaw(f1, e1, 0);
				sed.FitPowerlaw(f2, e2, 3);
				plotIStream << (f1/2.0) << " * (x*1000000)**" << e1 << " + " << (f2/2.0) << " * (x*1000000)**" << e2 << " with lines lw 1.0 title \"\"";
				
				if(sourceIndex != model.ComponentCount()-1)
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
				++sourceIndex;
			}
		}
	}
	
	if(!outputModel.empty())
		model.Save(outputModel.c_str());
	
	return 0;
}
