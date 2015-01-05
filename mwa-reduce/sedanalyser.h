#ifndef SED_ANALYSER_H
#define SED_ANALYSER_H

#include <fstream>
#include <memory>

#include "model.h"
#include "subbandpassband.h"
#include "progressbar.h"
#include "rmsynthesis.h"

class SEDAnalyser
{
public:
	SEDAnalyser() : _analysePassband(false) { }
	
	void Read(const char* filename)
	{
		_model.reset(new Model(filename));
	}
	
	void Set(const Model& model)
	{
		_model.reset(new Model(model));
	}
	
	long double CentralFrequency() const
	{
		return _model->begin()->begin()->SED().CentreFrequency();
	}
	
	const Model& GetModel() const { return *_model; }
	
	void Process();
	void ProcessRM();
	
	void OutputSIStatistics()
	{
		std::cout << "\nN_brightest_sources SI_avg SI_median SI_stddev SI_min SI_max lowestFlux\n";
		std::sort(_sources.rbegin(), _sources.rend());
		for(size_t n=10; n<_sources.size(); n*=2)
			outputSIStats(_sources, n);
		outputSIStats(_sources, _sources.size());
	}
	
	void SaveResults();
	
	double BestRMS();
	
	void SetAnalysePassband(bool analysePassband) { _analysePassband = analysePassband; }
	
private:
	std::unique_ptr<Model> _model;
	Model _brightestSourcesModel;
	bool _analysePassband;
	
	struct SourceInfo
	{
		size_t index, componentCount;
		long double plExponent, plFactor, pl2ndOrder;
		std::vector<double> terms;
		long double modalFlux, rms, rms2ndOrder, stokesVFrac, maxError, maxErrorFrequency;
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
		static bool hasLower2ndRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rms2ndOrder < rhs.rms2ndOrder; }

		std::string Description(const Model& model) const
		{
			std::ostringstream str;
			const ModelSource& ms = model.Source(index);
			const SpectralEnergyDistribution sed = ms.GetIntegratedSED();
			str << ms.Name() << ", " << modalFlux << " Jy, SI=" << plExponent << ", RMS=" << rms << " (" << rms2ndOrder << "), V%=" << round(10000.0*stokesVFrac)/100.0 << "%" << "(" << sed.AverageFlux(Polarization::StokesV) << "/" << sed.AverageFlux(Polarization::StokesI) << "), max E=" << maxError << " @ " << (maxErrorFrequency*1e-6) << ", RMpeak=" << rmPeakValue << " @ " << rmPeakPos << " / " << rmZeroValue << " @0, SB=" << subbandCorrelation << ", 2nd=" << pl2ndOrder;
			return str.str();
		}

	};
	
	std::vector<SourceInfo> _sources;
	
	double sourceRMS(const SpectralEnergyDistribution& sed, long double factor, long double exponent);
	
	double sourceRMS(const SpectralEnergyDistribution& sed, long double a, long double b, long double c);
	
	double sourceRMS(const SpectralEnergyDistribution& sed, const std::vector<double>& terms);
	
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
	
	void outputList(const std::vector<SourceInfo>& sortedList, const Model& model, const std::string& name, bool reverse);

	void outputSIStats(const std::vector<SourceInfo>& sortedList, size_t n);

};

#endif

