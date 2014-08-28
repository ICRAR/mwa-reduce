#include "banddata.h"
#include "model.h"
#include "imagecoordinates.h"
#include "areaset.h"
#include "parser/areaparser.h"
#include "loghistogram.h"
#include "beamevaluator.h"
#include "uvector.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <fstream>
#include <set>

class IsApparentLessBright
{
public:
	BeamEvaluator* _beam;
	double _frequency;
	
	bool operator()(const ModelSource& lhs, const ModelSource& rhs) const
	{
		return totalAppFlux(lhs) < totalAppFlux(rhs);
	}
	
	bool operator()(const ModelComponent& lhs, const ModelComponent& rhs) const
	{
		return totalAppFlux(lhs) < totalAppFlux(rhs);
	}
	
	double totalAppFlux(const ModelComponent& c) const
	{
		double stokesFlux[4];
		for(size_t p=0; p!=4; ++p)
			stokesFlux[p] = c.SED().FluxAtFrequency(_frequency, Polarization::IndexToStokes(p));
		std::complex<double> linFlux[4];
		Polarization::StokesToLinear(stokesFlux, linFlux);
		_beam->AbsToApparent(c.PosRA(), c.PosDec(), _frequency, linFlux);
		return (std::abs(linFlux[0]) + std::abs(linFlux[3])) * 0.5;
	}
	
	double totalAppFlux(const ModelSource& s) const
	{
		double totalFlux = 0.0;
		for(ModelSource::const_iterator compIter=s.begin(); compIter!=s.end(); ++compIter)
			totalFlux += totalAppFlux(*compIter);
		return totalFlux / s.ComponentCount();
	}
};

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		std::cout << "editmodel -- Interpolation, extrapolation, plotting and scaling of the \n"
		"spectral density function. Usage:\n"
		"\tsdf [-p] [-m <output model>] [-o] [-s <scale>] [-sp <peakflux A> <freq A> <peakflux B> <freq B>] [-sc <intflux A> <freq A> <intflux B> <freq B>] [-set0/1/2/3 <flux>] [-unpolarized] [-pl] [-t <threshold>] [-tapp <ms> <threshold>] [-r <new-nr-channels>] [-ravg <new-nr-channels>] [-near <ra> <dec> <dist>] [-combine-diff-meas] [<more models>..] [-collect <name>] [-rnd <n> <ra> <dec> <dist>] [-sort] [-appsort <ms>] [-without/only <areafile>] [-lognlogs <frequency> <bincount>] [-stats] [<model> [<more models...>]] [-delnans]\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, outputCsv = false, powerlaw = false, optimize = false, applyThreshold = false, applyAppThreshold = false, resample = false, resampleByAveraging = false, unpolarized = false, delNaNs = false;
	bool setPolarization[4] = {false, false, false, false};
	long double setPolFlux[4] = {0.0, 0.0, 0.0, 0.0};
	long double scale = 1.0, threshold = 0.0, appThreshold = 0.0, logNlogSFrequency = 0.0;
	long double scalePeakA = 1.0, scaleFreqA = 0.0, scalePeakB = 1.0, scaleFreqB = 0.0;
	size_t newChannelCount = 0, logNlogSBinCount = 0;
	std::string outputModel, collectName, appThresholdMS, appSortMS, csvFilename, plotTitle;
	bool nearFilter = false, scalePeak = false, scaleSource = false, doCollect = false, doSort = false, doAppSort = false;
	long double nearFilterRA = 0.0, nearFilterDec = 0.0, nearFilterDist = 0.0;
	enum { AddFluxes, AverageFluxes, DifferentFrequencies } combineStrategy = AddFluxes;
	bool excludeArea = false, logNlogS = false, sourceStats = false;
	std::unique_ptr<AreaSet> areaSet;
	size_t rndN = 0;
	double rndRA = 0.0, rndDec = 0.0, rndDist = 0.0;
	while(argi!=argc && argv[argi][0]=='-')
	{
		const std::string option(&argv[argi][1]);
		if(option == "p")
		{
			outputPlot = true;
		} else if(option == "ptitle") {
			++argi;
			outputPlot = true;
			plotTitle = argv[argi];
		} else if(option == "csv") {
			++argi;
			outputCsv = true;
			csvFilename = argv[argi];
		} else if(option == "near")
		{
			++argi;
			nearFilterRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			nearFilterDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			nearFilterDist = atof(argv[argi]) * (M_PI/180.0);
			nearFilter = true;
		} else if(option == "rnd")
		{
			++argi;
			rndN = atof(argv[argi]);
			++argi;
			rndRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			rndDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			rndDist = atof(argv[argi]) * (M_PI/180.0);
		} else if(option == "collect")
		{
			doCollect = true;
			++argi;
			collectName = argv[argi];
		} else if(option == "m")
		{
			++argi;
			outputModel = argv[argi];
		} else if(option == "s")
		{
			++argi;
			scale = atof(argv[argi]);
		} else if(option == "sp")
		{
			scalePeak = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(option == "sc")
		{
			scaleSource = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(option == "set0")
		{
			++argi;
			setPolarization[0] = true;
			setPolFlux[0] = atof(argv[argi]);
		} else if(option == "-set1")
		{
			++argi;
			setPolarization[1] = true;
			setPolFlux[1] = atof(argv[argi]);
		} else if(option == "set2")
		{
			++argi;
			setPolarization[2] = true;
			setPolFlux[2] = atof(argv[argi]);
		} else if(option == "set3")
		{
			++argi;
			setPolarization[3] = true;
			setPolFlux[3] = atof(argv[argi]);
		} else if(option == "pl")
		{
			powerlaw = true;
		} else if(option == "t")
		{
			++argi;
			threshold = atof(argv[argi]);
			applyThreshold = true;
		} else if(option == "tapp")
		{
			++argi;
			appThresholdMS = argv[argi];
			++argi;
			appThreshold = atof(argv[argi]);
			applyAppThreshold = true;
		} else if(option == "r")
		{
			++argi;
			newChannelCount = atoi(argv[argi]);
			resample = true;
		} else if(option == "ravg")
		{
			++argi;
			newChannelCount = atoi(argv[argi]);
			resampleByAveraging = true;
		} else if(option == "unpolarized")
		{
			unpolarized = true;
		} else if(option == "o")
		{
			optimize = true;
		} else if(option == "combine-diff-meas")
		{
			combineStrategy = DifferentFrequencies;
		} else if(option == "combine-avg-meas")
		{
			combineStrategy = AverageFluxes;
		} else if(option == "sort")
		{
			doSort = true;
		} else if(option == "appsort")
		{
			doAppSort = true;
			++argi;
			appSortMS = argv[argi];
		} else if(option == "without" || option == "only")
		{
			excludeArea = (option == "without");
			++argi;
			areaSet.reset(new AreaSet());
			AreaParser areaParser;
			areaParser.Parse(*areaSet, argv[argi]);
		} else if(option == "lognlogs")
		{
			++argi;
			logNlogSFrequency = atof(argv[argi]);
			++argi;
			logNlogSBinCount = atoi(argv[argi]);
			logNlogS = true;
		} else if(option == "stats")
		{
			sourceStats = true;
		} else if(option == "delnans")
		{
			delNaNs = true;
		} else {
			throw std::runtime_error(std::string("Unknown option given: -") + option);
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
						m->second.SetFluxDensityFromIndex(p, m->second.FluxDensityFromIndex(p) * fact);
				}
			}
		}
	}
	
	if(delNaNs)
	{
		for(Model::iterator s=model.begin(); s!=model.end(); ++s)
		{
			for(ModelSource::iterator c=s->begin(); c!=s->end(); ++c)
				c->SED().RemoveInvalidMeasurements();
		}
	}
	if(optimize)
		model.Optimize();
	if(doSort)
		model.Sort();
	if(doAppSort)
	{
		casa::MeasurementSet ms(appSortMS);
		BeamEvaluator beam(ms);
		BandData bandData(ms.spectralWindow());
		IsApparentLessBright compareFunc;
		compareFunc._beam = &beam;
		compareFunc._frequency = bandData.CentreFrequency();
		std::vector<ModelSource> temp(model.begin(), model.end());
		std::sort(temp.rbegin(), temp.rend(), compareFunc);
		
		Model tempModel;
		for(std::vector<ModelSource>::const_iterator s=temp.begin(); s!=temp.end(); ++s)
		{
			std::vector<ModelComponent> components(s->begin(), s->end());
			std::sort(components.rbegin(), components.rend(), compareFunc);
			ModelSource sortedSource(*s);
			sortedSource.ClearComponents();
			for(std::vector<ModelComponent>::const_iterator c=components.begin(); c!=components.end(); ++c)
				sortedSource.AddComponent(*c);
			tempModel.AddSource(sortedSource);
		}
		model = tempModel;
	}
	
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
			double flux = sourcePtr->TotalFlux(sed.LowestFrequency(), sed.HighestFrequency(), Polarization::StokesI);
			if(!std::isfinite(flux))
			{
				long double f, e;
				sed.FitPowerlaw(f, e, Polarization::StokesI);
				double lowFlux = f * std::pow(sed.LowestFrequency(), e);
				double hiFlux = f * std::pow(sed.HighestFrequency(), e);
				flux = (lowFlux + hiFlux) * 0.5;
			}

			if(flux >= threshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	if(applyAppThreshold)
	{
		Model temp;
		casa::MeasurementSet ms(appThresholdMS);
		BeamEvaluator beam(ms);
		BandData bandData(ms.spectralWindow());
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			double stokesFlux[4];
			for(size_t p=0; p!=4; ++p)
			{
				const SpectralEnergyDistribution& sed = sourcePtr->Peak().SED();
				stokesFlux[p] = sourcePtr->TotalFlux(bandData.LowestFrequency(), bandData.HighestFrequency(), Polarization::IndexToStokes(p));
				if(!std::isfinite(stokesFlux[p]))
				{
					long double f, e;
					sed.FitPowerlaw(f, e, Polarization::IndexToStokes(p));
					double lowFlux = f * std::pow(bandData.LowestFrequency(), e);
					double hiFlux = f * std::pow(bandData.HighestFrequency(), e);
					stokesFlux[p] = (lowFlux + hiFlux) * 0.5;
				}
			}
			std::complex<double> linFlux[4];
			Polarization::StokesToLinear(stokesFlux, linFlux);
			beam.AbsToApparent(sourcePtr->Peak().PosRA(), sourcePtr->Peak().PosDec(), bandData.CentreFrequency(), linFlux);
			double appStokesIFlux = (std::abs(linFlux[0]) + std::abs(linFlux[3])) * 0.5;
			if(appStokesIFlux >= appThreshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	if(rndN != 0)
	{
		for(size_t i=0; i!=rndN; ++i)
		{
			double ra, dec;
			do {
				ra = 2.0 * M_PI * (double) rand() / (double) RAND_MAX,
				dec = M_PI * ((double) rand() / (double) RAND_MAX - 0.5);
			} while(ImageCoordinates::AngularDistance(rndRA, rndDec, ra, dec) >= rndDist);
			ModelComponent component;
			component.SetPosRA(ra);
			component.SetPosDec(dec);
			component.SetSED(SpectralEnergyDistribution(1.0, 150000000.0));
			ModelSource source;
			source.AddComponent(component);
			std::ostringstream name;
			name << "rnd" << i;
			source.SetName(name.str());
			model.AddSource(source);
		}
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
		for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
		{
			const SpectralEnergyDistribution &sed = compPtr->SED();
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
					sed.FitPowerlaw(f, e, Polarization::IndexToStokes(p));
					m1.SetFluxDensityFromIndex(p, f * std::pow(startFreq, e));
					m2.SetFluxDensityFromIndex(p, f * std::pow(endFreq, e));
				}
				newSED.AddMeasurement(m1);
				newSED.AddMeasurement(m2);
				compPtr->SetSED(newSED);
			}
			else if(resample || resampleByAveraging)
			{
				SpectralEnergyDistribution newSED;
				const long double newBandSize = (endFreq - startFreq) / newChannelCount;
				for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
				{
					long double chStartFreq = startFreq + newBandSize*newChIndex;
					long double chEndFreq = startFreq + newBandSize*(newChIndex+1.0);
					
					Measurement m;
					m.SetFrequencyHz((chStartFreq+chEndFreq)*0.5);
					
					long double flux;
					if(newChannelCount >= sed.MeasurementCount())
						flux = sed.FluxAtFrequency((chStartFreq+chEndFreq)*0.5, Polarization::StokesI);
					else if(resampleByAveraging)
						flux = sed.AverageFlux(chStartFreq, chEndFreq, Polarization::StokesI);
					else 
						flux = sed.IntegratedFlux(chStartFreq, chEndFreq, Polarization::StokesI);
					m.SetFluxDensity(Polarization::StokesI, flux);
					
					PolarizationEnum pols[3] = { Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV };
					for(size_t p=0; p!=3; ++p)
					{
						if(newChannelCount >= sed.MeasurementCount())
							flux = sed.FluxAtFrequency((chStartFreq+chEndFreq)*0.5, pols[p]);
						else
							flux = sed.AverageFlux(chStartFreq, chEndFreq, pols[p]);
						m.SetFluxDensity(pols[p], flux);
					}
					
					newSED.AddMeasurement(m);
				}
				
				compPtr->SetSED(newSED);
			}
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
						m->second.SetFluxDensityFromIndex(p, setPolFlux[p]);
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
						m->second.SetFluxDensityFromIndex(p, m->second.FluxDensityFromIndex(p) * scale);
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
						oldFluxA = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqA, Polarization::IndexToStokes(p)),
						oldFluxB = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqB, Polarization::IndexToStokes(p));
					factorA[p] = oldFluxA==0.0 ? 0.0 : scalePeakA / oldFluxA;
					factorB[p] = oldFluxB==0.0 ? 0.0 : scalePeakB / oldFluxB;
				} else {
					long double
						oldFluxA = sourcePtr->TotalFlux(scaleFreqA, Polarization::IndexToStokes(p)),
						oldFluxB = sourcePtr->TotalFlux(scaleFreqB, Polarization::IndexToStokes(p));
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
					long double oldFluxA = compPtr->SED().FluxAtFrequency(scaleFreqA, Polarization::IndexToStokes(p));
					long double oldFluxB = compPtr->SED().FluxAtFrequency(scaleFreqB, Polarization::IndexToStokes(p));
					measA.SetFluxDensityFromIndex(p, oldFluxA*factorA[p]);
					measB.SetFluxDensityFromIndex(p, oldFluxB*factorB[p]);
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
	
	if(sourceStats)
	{
		const size_t n = model.SourceCount();
		double sum = 0.0, diffSq = 0.0;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double stokesI = source.TotalFlux(150000000.0, Polarization::StokesI);
			sum += stokesI;
			diffSq += (stokesI - 1.0) * (stokesI - 1.0);
		}
		double mean = sum / n, stdDevSum = 0.0;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double stokesI = source.TotalFlux(150000000.0, Polarization::StokesI);
			stdDevSum += (stokesI - mean) * (stokesI - mean);
		}
		std::cout <<
			"Source statistics:\n"
			" Mean: " << mean << "\n"
			" Standard deviation: " << sqrt(stdDevSum/n) << "\n"
			" Deviation from 1: " << sqrt(diffSq/n) << '\n';
	}
	
	if(logNlogS)
	{
		LogHistogram histogram;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double v = source.TotalFlux(logNlogSFrequency*1000000.0, Polarization::StokesI);
			histogram.Add(v);
		}
		long double max = histogram.MaxAmplitude(), min = histogram.MinPositiveAmplitude();
		long double step = powl(max/min, 1.0/logNlogSBinCount);
		std::ofstream plotDataStream("lognlogs-data.txt");
		ao::uvector<long double> vals(logNlogSBinCount);
		for(size_t i=0; i!=logNlogSBinCount; ++i)
		{
			long double
				intervalStart = min * powl(step, i),
				mid = min * powl(step, double(i) + 0.5),
				intervalEnd = min * powl(step, i+1),
				val = histogram.NormalizedCount(intervalStart, intervalEnd),
				valNormalized = log10(val * powl(mid, 2.5));
			vals[i] = valNormalized;
			plotDataStream << intervalStart << '\t' << intervalEnd << '\t' << intervalStart * sqrt(step) << '\t' << val << '\t' << valNormalized << '\n';
			intervalStart = intervalEnd;
		}
		size_t avgCount = 0;
		long double avg = 0.0;
		for(size_t i=logNlogSBinCount/5*3; i<logNlogSBinCount; ++i)
		{
			if(std::isfinite(vals[i]))
			{
				avg += vals[i];
				++avgCount;
			}
		}
		avg /= avgCount;
		
		std::ofstream plotStream("lognlogs.plt");
		plotStream <<
			"set terminal postscript enhanced color\n"
			"set logscale x\n"
			"set xrange [" << min << ":" << max << "]\n"
			"set yrange [:1.2]\n"
			"set output \"lognlogs.ps\"\n"
			"set key bottom right\n"
			"set xlabel \"Flux (Jy)\"\n"
			"set ylabel \"log N/S^{-2.5}/avg\"\n"
			"plot \"lognlogs-data.txt\" using 2:(column(5)/" << avg << ") with lines lw 2.0 title \"Measured\", \\\n"
			"1 with lines lw 2.0 lt 3 title \"Euclidean\"\n";
	}
	
	if(outputCsv)
	{
		std::ofstream csvFile(csvFilename);
		csvFile << "freq,i,q,u,v\n";
		if(model.SourceCount() != 0)
		{
			if(model.SourceCount() != 1)
				std::cout << "Warning: multiple sources in model, but will only output a csv file for the first source.\n";
			const ModelSource& source = *model.begin();
			if(source.ComponentCount() != 1)
				std::cout << "Warning: first source has multiple components; only outputting first component.\n";
			const SpectralEnergyDistribution &sed = source.begin()->SED();

			std::vector<Measurement> measurements;
			sed.GetMeasurements(measurements);
			
			for(SpectralEnergyDistribution::const_iterator iter=sed.begin(); iter!=sed.end(); ++iter)
			{
				const Measurement& m = iter->second;
				long double
					i = m.FluxDensity(Polarization::StokesI), q = m.FluxDensity(Polarization::StokesQ),
					u = m.FluxDensity(Polarization::StokesU), v = m.FluxDensity(Polarization::StokesV);
				if(std::isfinite(i) && std::isfinite(q) && std::isfinite(u) && std::isfinite(v))
				{
					csvFile
						<< m.FrequencyHz()*1e-6 << ','
						<< i << ',' << q << ',' << u << ',' << v << '\n';
				}
			}
		}
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
			"#set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n";
		if(!plotTitle.empty())
			plotStream << "set title \"" << plotTitle << "\"\n";
		plotStream <<
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
			"set ylabel \"Flux (Jy)\"\n";
		if(!plotTitle.empty())
			plotIStream << "set title \"" << plotTitle << "\"\n";
		plotIStream <<
			"plot \\\n";
		std::ofstream plot10IStream("spectrum10-I.plt");
		plot10IStream <<
			"set terminal postscript enhanced color\n"
			"set logscale y\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"spectrum10-I.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"Frequency (MHz)\"\n"
			"set ylabel \"Flux (Jy)\"\n";
		if(!plotTitle.empty())
			plot10IStream << "set title \"" << plotTitle << "\"\n";
		plot10IStream <<
			"plot \\\n";

		size_t compIndex = 0, sourceIndex = 0;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				std::ostringstream dataStreamName;
				dataStreamName << "spectrum" << compIndex << ".txt";
				std::ofstream dataStream(dataStreamName.str().c_str());
				if(compIndex == 0)
				{
					plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with points lw 1.0 title \"I\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:3 with points lw 1.0 title \"Q\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:4 with points lw 1.0 title \"U\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:5 with points lw 1.0 title \"V\"";
				}
				else {
					plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with points lw 1.0 title \"\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:3 with points lw 1.0 title \"\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:4 with points lw 1.0 title \"\",\\\n";
					plotStream << "\"" << dataStreamName.str() << "\" using 1:5 with points lw 1.0 title \"\"";
				}
				
				plotIStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 1.0 title \"\",\\\n";
				if(sourceIndex < 10 && compPtr == sourcePtr->begin())
					plot10IStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 1.0 title \"" << sourcePtr->Name() << "\",\\\n";
				
				const SpectralEnergyDistribution &sed = compPtr->SED();
				/*long double e, f;
				sed.FitPowerlaw(f, e, Polarization::StokesI);
				plotIStream << (f/2.0) << " * (x*1000000)**" << e << " with lines lw 1.0 title \"\"";*/
				
				if(compIndex != model.ComponentCount()-1)
				{
					plotStream << ",\\";
					//plotIStream << ",\\";
				}
				plotStream << "\n";
				//plotIStream << "\n";
				
				std::vector<Measurement> measurements;
				sed.GetMeasurements(measurements);
				
				for(std::vector<Measurement>::const_iterator i=measurements.begin(); i!=measurements.end(); ++i)
				{
					dataStream
						<< i->FrequencyHz()/1000000.0 << '\t'
						<< i->FluxDensity(Polarization::StokesI) << '\t'
						<< i->FluxDensity(Polarization::StokesQ) << '\t'
						<< i->FluxDensity(Polarization::StokesU) << '\t'
						<< i->FluxDensity(Polarization::StokesV) << '\n';
				}
				++compIndex;
			}
			++sourceIndex;
		}
	}
	
	if(!outputModel.empty()) {
		std::cout << "Writing model with " << model.SourceCount() << " sources.\n";
		model.Save(outputModel.c_str());
	}
	
	return 0;
}
