#include "sedanalyser.h"

void SEDAnalyser::Process()
{
	const Model& model = *_model;
	long double centreFrequency = CentralFrequency();
	ProgressBar progress("Analysing SEDs");
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const ModelSource& source = model.Source(i);
		SourceInfo newSource;
		const SpectralEnergyDistribution sed = source.GetIntegratedSED();
		sed.FitPowerlaw(newSource.plFactor, newSource.plExponent, Polarization::StokesI);
		newSource.index = i;
		newSource.modalFlux = newSource.plFactor * powl(centreFrequency, newSource.plExponent);
		newSource.componentCount = source.ComponentCount();
		newSource.rms = sourceRMS(sed, newSource.plFactor, newSource.plExponent);
		newSource.stokesVFrac = std::fabs(sed.AverageFlux(Polarization::StokesV) / sed.AverageFlux(Polarization::StokesI));
		newSource.maxError = 0.0;
		for(SpectralEnergyDistribution::const_iterator m=sed.begin(); m!=sed.end(); ++m)
		{
			// There seems to be some consistent high values in the last two channels... Systematics?
			// Ignore for now
			if(m->first > 197.61*1e6)
				break;
			double err = (m->second.FluxDensity(Polarization::StokesI) - newSource.plFactor * powl(m->first, newSource.plExponent)) / newSource.rms;
			if(fabs(err) > fabs(newSource.maxError))
			{
				newSource.maxError = err;
				newSource.maxErrorFrequency = m->first;
			}
		}
		long double aTemp, bTemp;
		sed.FitPowerlaw2ndOrder(aTemp, bTemp, newSource.pl2ndOrder, Polarization::StokesI);
		newSource.rms2ndOrder = sourceRMS(sed, aTemp, bTemp, newSource.pl2ndOrder);
		//std::cout << "First fit=" << newSource.plExponent << ", 2nd=" << aTemp << ", c=" << newSource.pl2ndOrder << '\n';
		newSource.subbandCorrelation = getSubbandShapeCorrelation(sed);
		if(source.ComponentCount() == 1) // skip complex sources
			_sources.push_back(newSource);
		
		progress.SetProgress(i+1, model.SourceCount());
	}
}

void SEDAnalyser::ProcessRM()
{
	const Model& model = *_model;
	ProgressBar progress("Doing RM synthesis");
	for(size_t i=0; i!=_sources.size(); ++i)
	{
		SourceInfo& newSource = _sources[i];
		const ModelSource& source = model.Source(newSource.index);
		const SpectralEnergyDistribution sed = source.GetIntegratedSED();
		RMSynthesis rmSynthesis(sed);
		rmSynthesis.Synthesize();
		rmSynthesis.PeakWithNonZeroRM(5, newSource.rmPeakPos, newSource.rmPeakValue);
		double tmp;
		rmSynthesis.PeakWithSmallRM(5, tmp, newSource.rmZeroValue);
		progress.SetProgress(i+1, model.SourceCount());
	}
}

void SEDAnalyser::SaveResults()
{

	{
		std::cout << "Sorting... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend());
		std::cout << "DONE\n";
		
		std::cout << "Name\tcomplex?\tS" << CentralFrequency()*1e-6 << "\tSI\tRMS\n";
		
		std::vector<std::pair<size_t,double>> siAverages;
		double spectralIndexSum = 0.0;
		size_t siAvgCount = 0;
		
		const Model& model = *_model;
		_brightestSourcesModel = Model();
		for(std::vector<SourceInfo>::const_iterator i=_sources.begin(); i!=_sources.end(); ++i)
		{
			const ModelSource& ms = model.Source(i->index);
			std::cout << ms.Name() << '\t' << i->modalFlux << '\t' << i->plExponent << '\t' << i->rms << '\n';
			spectralIndexSum += i->plExponent;
			++siAvgCount;
			if(siAvgCount < 20) _brightestSourcesModel.AddSource(ms);
		}
	
		std::sort(_sources.rbegin(), _sources.rend());
		_brightestSourcesModel.Save("brightest-sources-model.txt");
		
		std::cout << "\nSorting on RMS... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerRMS);
		std::cout << "DONE\n";
		std::cout << "High RMS sources:\n";
		outputList(_sources, *_model, "highrms", false);
		std::cout << "\nLow RMS sources:\n";
		outputList(_sources, *_model, "lowrms", true);
		
		std::cout << "\nSorting on spectral index... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerSI);
		std::cout << "DONE\n";
		std::cout << "High SI sources:\n";
		outputList(_sources, *_model, "highsi", false);
		std::cout << "\nLow SI sources:\n";
		outputList(_sources, *_model, "lowsi", true);
		
		std::cout << "\nSorting on Stokes V frac... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerVFrac);
		std::cout << "DONE\n";
		std::cout << "High V-fractional sources:\n";
		outputList(_sources, *_model, "highvfrac", false);

		std::cout << "\nSorting on max error... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerMaxError);
		std::cout << "DONE\n";
		std::cout << "Max errors:\n";
		outputList(_sources, *_model, "maxerror", false);
		
		std::cout << "\nSorting on RM value... " << std::flush;
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerRMValue);
		std::cout << "DONE\n";
		std::cout << "Relatively high RM values:\n";
		outputList(_sources, *_model, "highrm", false);
		
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerSubbandCorrelation);
		std::cout << "High sub-band correlation:\n";
		outputList(_sources, *_model, "highsubbandcorrelation", false);
		std::cout << "Low sub-band correlation:\n";
		outputList(_sources, *_model, "lowsubbandcorrelation", true);
		
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLower2ndOrder);
		std::cout << "High 2nd order value:\n";
		outputList(_sources, *_model, "high2ndorder", false);
		std::cout << "Low 2nd order value:\n";
		outputList(_sources, *_model, "low2ndorder", true);
		
		std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLower2ndRMSDiff);
		std::cout << "High 2nd order RMS diff:\n";
		outputList(_sources, *_model, "high2ndrmsdiff", false);
		std::cout << "Low 2nd order RMS diff:\n";
		outputList(_sources, *_model, "low2ndrmsdiff", true);
	}
}

void SEDAnalyser::outputList(const std::vector< SEDAnalyser::SourceInfo >& sortedList, const Model& model, const std::string& name, bool reverse)
{
	{
		Model outputModel;
		for(size_t i=0; i<std::min(size_t(20), sortedList.size()); ++i)
		{
			size_t index = reverse ? (sortedList.size()-i-1) : i;
			const SourceInfo& source = sortedList[index];
			outputModel.AddSource(model.Source(source.index));
			std::cout << source.Description(model) << '\n';
		}
		outputModel.Save(name + "-sources-model.txt");
	}
}

void SEDAnalyser::outputSIStats(const std::vector< SEDAnalyser::SourceInfo >& sortedList, size_t n)
{

	{
		size_t actualN = n > sortedList.size() ? sortedList.size() : n;
		std::vector<long double> sis(actualN);
		long double sum = 0.0;
		for(size_t i=0; i!=actualN; ++i)
		{
			const SourceInfo& source = sortedList[i];
			sis[i] = source.plExponent;
			sum += source.plExponent;
		}
		long double lowestFlux = sortedList[actualN-1].modalFlux;
		long double average = sum / actualN;
		long double squaredDistSum = 0.0;
		for(size_t i=0; i!=actualN; ++i)
		{
			long double dist = sis[i]-average;
			squaredDistSum += dist*dist;
		}
		long double stddev = sqrtl(squaredDistSum/actualN);
		
		std::sort(sis.begin(), sis.end());
		long double median;
		if(actualN%2 == 0)
			median = (sis[actualN/2-1] + sis[actualN/2]) * 0.5;
		else
			median = sis[actualN/2];
		long double highestSI = sis.back(), lowestSI = sis.front();
		
		std::cout << actualN << ' ' << average << ' ' << median << ' ' << stddev << ' ' << lowestSI << ' ' << highestSI << ' ' << lowestFlux << ' ' << '\n';
		
		std::vector<size_t> histData(101, 0);
		for(size_t i=0; i!=actualN; ++i)
		{
			long double bindex = ((long double) histData.size())*(sis[i]-lowestSI) / (highestSI-lowestSI);
			histData[size_t(round(bindex))]++;
		}
		
		std::ostringstream histFilename;
		histFilename << "histogram-n" << actualN << ".txt";
		std::ofstream histFile(histFilename.str());
		for(size_t i=0; i!=histData.size(); ++i)
		{
			long double binCentre = ((long double) (i))/histData.size()*(highestSI-lowestSI)+lowestSI;
			histFile << i << '\t' << binCentre << '\t' << histData[i] << '\n';
		}
	}
}

double SEDAnalyser::BestRMS()
{
	std::sort(_sources.begin(), _sources.end(), &SourceInfo::hasLower2ndRMS);
	const SourceInfo& source = _sources.front();
	return source.rms2ndOrder;
}
