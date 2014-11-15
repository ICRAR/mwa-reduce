#include <iostream>
#include <fstream>

#include "model.h"
#include "rmsynthesis.h"
#include "progressbar.h"
#include "subbandpassband.h"

struct SourceInfo
{
	size_t index, componentCount;
	long double plExponent, plFactor, pl2ndOrder;
	long double modalFlux, rms, rms2ndOrder, stokesVFrac, maxError;
	double rmPeakValue, rmPeakPos, rmZeroValue, subbandCorrelation;
	
	double RMRelValue() const { return rmPeakPos / rmZeroValue; }
	
	bool operator<(const SourceInfo& rhs) const { return modalFlux < rhs.modalFlux; }
	static bool hasLowerRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rms < rhs.rms; }
	static bool hasLowerSI(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.plExponent < rhs.plExponent; }
	static bool hasLowerVFrac(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.stokesVFrac < rhs.stokesVFrac; }
	static bool hasLowerMaxError(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.maxError < rhs.maxError; }
	static bool hasLowerRMValue(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rmPeakValue < rhs.rmPeakValue; }
	static bool hasLowerRMRelValue(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.RMRelValue() < rhs.RMRelValue(); }
	static bool hasLowerSubbandCorrelation(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.subbandCorrelation < rhs.subbandCorrelation; }
	static bool hasLower2ndOrder(const SourceInfo& lhs, const SourceInfo& rhs) { return std::fabs(lhs.pl2ndOrder) < std::fabs(rhs.pl2ndOrder); }
	static bool hasLower2ndRMSDiff(const SourceInfo& lhs, const SourceInfo& rhs) { return std::fabs(lhs.rms - lhs.rms2ndOrder) < std::fabs(rhs.rms - rhs.rms2ndOrder); }
};

double SourceRMS(const SpectralEnergyDistribution& sed, long double factor, long double exponent)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(SpectralEnergyDistribution::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(Polarization::StokesI);
		long double powerLawValue = factor * powl(i->second.FrequencyHz(), exponent);
		if(std::isfinite(flux))
		{
			long double diff = flux - powerLawValue;
			++count;
			RMSsum += diff * diff;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

double SourceRMS(const SpectralEnergyDistribution& sed, long double a, long double b, long double c)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(SpectralEnergyDistribution::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(Polarization::StokesI);
		long double x = i->second.FrequencyHz();
		long double powerLawValue = powl(b*x + c*x*x, a);
		if(std::isfinite(flux))
		{
			long double diff = flux - powerLawValue;
			++count;
			RMSsum += diff * diff;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

std::string SourceDesc(const SourceInfo& source, const Model& model)
{
	std::ostringstream str;
	const ModelSource& ms = model.Source(source.index);
	const SpectralEnergyDistribution sed = ms.GetIntegratedSED();
	str << ms.Name() << ", " << source.modalFlux << " Jy, SI=" << source.plExponent << ", RMS=" << source.rms << " (" << source.rms2ndOrder << "), V%=" << round(10000.0*source.stokesVFrac)/100.0 << "%" << "(" << sed.AverageFlux(Polarization::StokesV) << "/" << sed.AverageFlux(Polarization::StokesI) << "), max E=" << source.maxError << ", RMpeak=" << source.rmPeakValue << " @ " << source.rmPeakPos << " / " << source.rmZeroValue << " @0, SB=" << source.subbandCorrelation << ", 2nd=" << source.pl2ndOrder;
	return str.str();
}

void outputList(const std::vector<SourceInfo>& sortedList, const Model& model, const std::string& name, bool reverse)
{
	Model outputModel;
	for(size_t i=0; i<std::min(size_t(20), sortedList.size()); ++i)
	{
		size_t index = reverse ? (sortedList.size()-i-1) : i;
		const SourceInfo& source = sortedList[index];
		outputModel.AddSource(model.Source(source.index));
		std::cout << SourceDesc(source, model) << '\n';
	}
	outputModel.Save(name + "-sources-model.txt");
}


void makeMeanZero(std::vector<double>& values)
{
	size_t count = 0;
	double avg = 0.0;
	for(size_t i=0; i!=values.size(); ++i)
	{
		if(std::isfinite(values[i]))
		{
			avg += values[i];
			++count;
		}
	}
	avg /= count;
	for(size_t i=0; i!=values.size(); ++i)
		values[i] -= avg;
}

double getSubbandShapeCorrelation(std::vector<double>& measurements, const std::vector<double>& normalizedPassband)
{
	makeMeanZero(measurements);
	if(measurements.size() != normalizedPassband.size())
		throw std::runtime_error("Number of channels in passband and number of channels in selected part of measurements do not match");
	double match = 0.0;
	size_t count = 0;
	for(size_t i=0; i!=normalizedPassband.size(); ++i)
	{
		if(std::isfinite(normalizedPassband[i]) && std::isfinite(measurements[i]))
		{
			++count;
			match += measurements[i] * normalizedPassband[i];
		}
	}
	return match / count;
}

double getSubbandShapeCorrelation(const SpectralEnergyDistribution& sed)
{
	std::vector<double> passband;
	size_t channelsInPassband = 32;
	SubbandPassband::GetSubbandPassband(passband, channelsInPassband);
	makeMeanZero(passband);
	
	double bandStart = sed.begin()->second.FrequencyHz();
	std::vector<double> measurements;
	double correlationSum = 0.0;
	size_t correlationCount = 1;
	for(SpectralEnergyDistribution::const_iterator m = sed.begin(); m!=sed.end(); ++m)
	{
		const Measurement& meas = m->second;
		if(meas.FrequencyHz() - bandStart > 1280000.0 - 640000.0 / (measurements.size()+1))
		{
			correlationSum += std::fabs(getSubbandShapeCorrelation(measurements, passband));
			++correlationCount;
			measurements.clear();
			bandStart = meas.FrequencyHz();
		}
		measurements.push_back(meas.FluxDensity(Polarization::StokesI));
	}
	correlationSum += std::fabs(getSubbandShapeCorrelation(measurements, passband));
	return correlationSum / correlationCount;
}

void outputSIStats(const std::vector<SourceInfo>& sortedList, size_t n)
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

int main(int argc, char* argv[])
{
	bool withRM = false;
	
	std::cout << "Reading model... " << std::flush;
	const char* filename = argv[1];
	Model model(filename);
	std::cout << "DONE\n";
	
	std::cout << "Determining centre frequency... " << std::flush;
	long double centreFrequency = model.begin()->begin()->SED().CentreFrequency();
	std::cout << centreFrequency*1e-6 << " MHz.\n";
	
	std::cout << "Fitting power laws... " << std::flush;
	std::vector<SourceInfo> sources;
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const ModelSource& source = model.Source(i);
		SourceInfo newSource;
		const SpectralEnergyDistribution sed = source.GetIntegratedSED();
		sed.FitPowerlaw(newSource.plFactor, newSource.plExponent, Polarization::StokesI);
		newSource.index = i;
		newSource.modalFlux = newSource.plFactor * powl(centreFrequency, newSource.plExponent);
		newSource.componentCount = source.ComponentCount();
		newSource.rms = SourceRMS(sed, newSource.plFactor, newSource.plExponent);
		newSource.stokesVFrac = std::fabs(sed.AverageFlux(Polarization::StokesV) / sed.AverageFlux(Polarization::StokesI));
		newSource.maxError = 0.0;
		for(SpectralEnergyDistribution::const_iterator m=sed.begin(); m!=sed.end(); ++m)
		{
			double err = (m->second.FluxDensity(Polarization::StokesI) - newSource.plFactor * powl(m->first, newSource.plExponent)) / newSource.rms;
			if(err > newSource.maxError) newSource.maxError = err;
		}
		long double aTemp, bTemp;
		sed.FitPowerlaw2ndOrder(aTemp, bTemp, newSource.pl2ndOrder, Polarization::StokesI);
		newSource.rms2ndOrder = SourceRMS(sed, aTemp, bTemp, newSource.pl2ndOrder);
		std::cout << "First fit=" << newSource.plExponent << ", 2nd=" << aTemp << ", c=" << newSource.pl2ndOrder << '\n';
		newSource.subbandCorrelation = getSubbandShapeCorrelation(sed);
		if(source.ComponentCount() == 1) // skip complex sources
			sources.push_back(newSource);
	}
	std::cout << "DONE\n";
	
	if(withRM)
	{
		ProgressBar progress("Doing RM synthesis");
		for(size_t i=0; i!=sources.size(); ++i)
		{
			SourceInfo& newSource = sources[i];
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
	
	std::cout << "Sorting... " << std::flush;
	std::sort(sources.rbegin(), sources.rend());
	std::cout << "DONE\n";
	
	std::cout << "Name\tcomplex?\tS" << centreFrequency*1e-6 << "\tSI\tRMS\n";
	
	std::vector<std::pair<size_t,double>> siAverages;
	double spectralIndexSum = 0.0;
	size_t siAvgCount = 0;
	Model brightestSourcesModel;
	for(std::vector<SourceInfo>::const_iterator i=sources.begin(); i!=sources.end(); ++i)
	{
		const ModelSource& ms = model.Source(i->index);
		std::cout << ms.Name() << '\t' << i->modalFlux << '\t' << i->plExponent << '\t' << i->rms << '\n';
		spectralIndexSum += i->plExponent;
		++siAvgCount;
		if(siAvgCount < 20) brightestSourcesModel.AddSource(ms);
	}
	brightestSourcesModel.Save("brightest-sources-model.txt");
	
	std::cout << "\nN_brightest_sources SI_avg SI_median SI_stddev SI_min SI_max lowestFlux\n";
	for(size_t n=10; n<sources.size(); n*=2)
		outputSIStats(sources, n);
	outputSIStats(sources, sources.size());
	
	std::cout << "\nSorting on RMS... " << std::flush;
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerRMS);
	std::cout << "DONE\n";
	std::cout << "High RMS sources:\n";
	outputList(sources, model, "highrms", false);
	std::cout << "\nLow RMS sources:\n";
	outputList(sources, model, "lowrms", true);
	
	std::cout << "\nSorting on spectral index... " << std::flush;
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerSI);
	std::cout << "DONE\n";
	std::cout << "High SI sources:\n";
	outputList(sources, model, "highsi", false);
	std::cout << "\nLow SI sources:\n";
	outputList(sources, model, "lowsi", true);
	
	std::cout << "\nSorting on Stokes V frac... " << std::flush;
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerVFrac);
	std::cout << "DONE\n";
	std::cout << "High V-fractional sources:\n";
	outputList(sources, model, "highvfrac", false);

	std::cout << "\nSorting on max error... " << std::flush;
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerMaxError);
	std::cout << "DONE\n";
	std::cout << "Max errors:\n";
	outputList(sources, model, "maxerror", false);
	
	std::cout << "\nSorting on RM value... " << std::flush;
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerRMValue);
	std::cout << "DONE\n";
	std::cout << "Relatively high RM values:\n";
	outputList(sources, model, "highrm", false);
	
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLowerSubbandCorrelation);
	std::cout << "High sub-band correlation:\n";
	outputList(sources, model, "highsubbandcorrelation", false);
	std::cout << "Low sub-band correlation:\n";
	outputList(sources, model, "lowsubbandcorrelation", true);
	
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLower2ndOrder);
	std::cout << "High 2nd order value:\n";
	outputList(sources, model, "high2ndorder", false);
	std::cout << "Low 2nd order value:\n";
	outputList(sources, model, "low2ndorder", true);
	
	std::sort(sources.rbegin(), sources.rend(), &SourceInfo::hasLower2ndRMSDiff);
	std::cout << "High 2nd order RMS diff:\n";
	outputList(sources, model, "high2ndrmsdiff", false);
	std::cout << "Low 2nd order RMS diff:\n";
	outputList(sources, model, "low2ndrmsdiff", true);
}
