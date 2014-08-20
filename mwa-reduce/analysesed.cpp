#include <iostream>

#include "model.h"

struct SourceInfo
{
	size_t index, componentCount;
	long double plExponent, plFactor;
	long double modalFlux, rms, stokesVFrac;
	
	bool operator<(const SourceInfo& rhs) const { return modalFlux < rhs.modalFlux; }
	static bool hasLowerRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rms < rhs.rms; }
	static bool hasLowerSI(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.plExponent < rhs.plExponent; }
	static bool hasLowerVFrac(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.stokesVFrac < rhs.stokesVFrac; }
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

std::string SourceDesc(const SourceInfo& source, const Model& model)
{
	std::ostringstream str;
	const ModelSource& ms = model.Source(source.index);
	const SpectralEnergyDistribution sed = ms.GetIntegratedSED();
	str << ms.Name() << ", " << source.modalFlux << " Jy, SI=" << source.plExponent << ", RMS=" << source.rms << ", V%=" << round(10000.0*source.stokesVFrac)/100.0 << "%" << "(" << sed.AverageFlux(Polarization::StokesV) << "/" << sed.AverageFlux(Polarization::StokesI) << ")";
	return str.str();
}

void outputList(const std::vector<SourceInfo>& sortedList, const Model& model, const std::string& name, bool reverse)
{
	Model outputModel;
	for(size_t i=0; i<std::min(size_t(50), sortedList.size()); ++i)
	{
		size_t index = reverse ? (sortedList.size()-i-1) : i;
		const SourceInfo& source = sortedList[index];
		outputModel.AddSource(model.Source(source.index));
		std::cout << SourceDesc(source, model) << '\n';
	}
	outputModel.Save(name + "-sources-model.txt");
}

int main(int argc, char* argv[])
{
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
		if(source.ComponentCount() == 1) // skip complex sources
			sources.push_back(newSource);
	}
	std::cout << "DONE\n";
	
	std::cout << "Sorting... " << std::flush;
	std::sort(sources.rbegin(), sources.rend());
	std::cout << "DONE\n";
	
	std::cout << "Name\tcomplex?\tS" << centreFrequency*1e-6 << "\tSI\tRMS\n";
	
	std::vector<std::pair<size_t,double>> siAverages;
	double spectralIndexSum = 0.0;
	size_t siAvgCount = 0, nextSICount = 10;
	Model brightestSourcesModel;
	for(std::vector<SourceInfo>::const_iterator i=sources.begin(); i!=sources.end(); ++i)
	{
		const ModelSource& ms = model.Source(i->index);
		std::cout << ms.Name() << '\t' << i->modalFlux << '\t' << i->plExponent << '\t' << i->rms << '\n';
		spectralIndexSum += i->plExponent;
		++siAvgCount;
		if(siAvgCount == nextSICount)
		{
			siAverages.push_back(std::make_pair(siAvgCount, spectralIndexSum/(long double)(siAvgCount)));
			nextSICount *= 2;
		}
		if(siAvgCount < 20) brightestSourcesModel.AddSource(ms);
	}
	brightestSourcesModel.Save("brightest-sources-model.txt");
	
	std::cout << "\nN_brightest_sources\tavg_SI\n";
	for(std::vector<std::pair<size_t,double>>::const_iterator i=siAverages.begin(); i!=siAverages.end(); ++i)
		std::cout << i->first << '\t' << i->second << '\n';
	std::cout << siAvgCount << '\t' << (spectralIndexSum/(long double)(siAvgCount)) << '\n';
	
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
}
