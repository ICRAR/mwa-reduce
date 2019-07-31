#ifndef SED_ANALYSER_H
#define SED_ANALYSER_H

#include <fstream>
#include <memory>
#include <random>

#include "../../model/model.h"
#include "../../progressbar.h"
#include "../../rmsynthesis.h"
#include "../../subbandpassband.h"

class SEDAnalyser
{
public:
	SEDAnalyser() : _analysePassband(false), _edgeChannels(3), _flaggedEndChannels(0), _siBinCount(4), _powerSpectrumAnalysis(false), _sourceCountAnalysis(false), _useLombPeriodogram(false), _noMissing(false), _psInTemperature(false), _referenceFrequency(1), _nTermsInPS(3) { }
	
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
		if(_model==0 || _model->Empty())
			return std::numeric_limits<long double>::quiet_NaN();
		else
			return _model->begin()->begin()->MSED().CentreFrequency();
	}
	
	const Model& GetModel() const { return *_model; }
	
	void Process();
	void ProcessRM();
	
	void OutputSIStatistics()
	{
		outputGeneralSIStats();
		std::cout << "\nN_brightest_sources SI_avg SI_median SI_stddev SI_min SI_max lowestFlux\n";
		std::sort(_sources.rbegin(), _sources.rend());
		for(size_t n=10; n<_sources.size(); n*=2)
			outputSIStats(_sources, n);
		outputSIStats(_sources, _sources.size());
	}
	
	void OutputSIBrightnessStatistics();
	
	void SaveResults();
	
	double BestRMS();
	
	void SetAnalysePassband(bool analysePassband) { _analysePassband = analysePassband; }
	
	void SetEdgeChannelCount(size_t edgeChannels) { _edgeChannels = edgeChannels; }
	
	void SetFlaggedEndChannelCount(size_t flaggedEndChannels) { _flaggedEndChannels = flaggedEndChannels; }
	
	void SetPowerSpectrumAnalysis(bool psAnalysis) { _powerSpectrumAnalysis = psAnalysis; }
	
	void SetSourceCountAnalysis(bool sourceCountAnalysis) { _sourceCountAnalysis = sourceCountAnalysis; }
	
	void SetLombPeriodogramAnalysis(bool lombAnalysis) { _useLombPeriodogram = lombAnalysis; }
	
	void SetNoMissing(bool noMissing) { _noMissing = noMissing; }
	
	void SetPSInTemperature(bool useTemperature) { _psInTemperature = useTemperature; }
	
	void SetSIBinCount(size_t binCount) { _siBinCount = binCount; }
	
	void SetReferenceFrequency(double refFreq) { _referenceFrequency = refFreq; }
	
	void SetTermsInPS(size_t terms) { _nTermsInPS = terms; }
	
	void MakePSPlot();
	void MakeAveragedPSPlot();
	void MakeNoisePSPlot();
private:
	std::unique_ptr<Model> _model;
	Model _brightestSourcesModel;
	bool _analysePassband;
	size_t _edgeChannels, _flaggedEndChannels, _siBinCount;
	bool _powerSpectrumAnalysis, _sourceCountAnalysis, _useLombPeriodogram, _noMissing, _psInTemperature;
	double _referenceFrequency;
	size_t _nTermsInPS;
	
	struct SourceInfo
	{
		size_t index, componentCount;
		long double plExponent, plFactor, pl2ndOrder;
		ao::uvector<double> terms;
		long double modalFlux, rms, qrms, urms, vrms, rms2ndOrder, diffRMS, stokesVFrac, maxError, maxAbsError, minAbsError, maxErrorFrequency, maxStep, distance;
		double rmPeakValue, rmPeakPos, rmZeroValue, subbandCorrelation;
		
		double RMRelValue() const { return rmPeakPos / rmZeroValue; }
		
		bool operator<(const SourceInfo& rhs) const { return modalFlux < rhs.modalFlux; }
		static bool hasLowerRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rms < rhs.rms; }
		static bool hasLowerSI(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.plExponent < rhs.plExponent; }
		static bool hasLowerVFrac(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.stokesVFrac < rhs.stokesVFrac; }
		static bool hasLowerMaxError(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.maxError < rhs.maxError; }
		static bool hasLowerMaxAbsError(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.maxAbsError < rhs.maxAbsError; }
		static bool hasLowerMinAbsError(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.minAbsError < rhs.minAbsError; }
		static bool hasLowerRMValue(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rmPeakValue < rhs.rmPeakValue; }
		static bool hasLowerRMRelValue(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.RMRelValue() < rhs.RMRelValue(); }
		static bool hasLowerSubbandCorrelation(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.subbandCorrelation < rhs.subbandCorrelation; }
		static bool hasLower2ndOrder(const SourceInfo& lhs, const SourceInfo& rhs) { return std::fabs(lhs.pl2ndOrder) < std::fabs(rhs.pl2ndOrder); }
		static bool hasLower2ndRMSDiff(const SourceInfo& lhs, const SourceInfo& rhs) { return std::fabs(lhs.rms - lhs.rms2ndOrder) < std::fabs(rhs.rms - rhs.rms2ndOrder); }
		static bool hasLower2ndRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.rms2ndOrder < rhs.rms2ndOrder; }
		static bool hasLowerDiffRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.diffRMS < rhs.diffRMS; }
		static bool hasLowerQRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.qrms < rhs.qrms; }
		static bool hasLowerURMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.urms < rhs.urms; }
		static bool hasLowerVRMS(const SourceInfo& lhs, const SourceInfo& rhs) { return lhs.vrms < rhs.vrms; }

		std::string Description(const Model& model) const
		{
			std::ostringstream str;
			const ModelSource& ms = model.Source(index);
			const MeasuredSED sed = ms.GetIntegratedMSED();
			str << ms.Name() << ", "
				<< RaDecCoord::RAToString(ms.MeanRA()) << ' '
				<< RaDecCoord::DecToString(ms.MeanDec()) << ", "
				<< modalFlux << " Jy, SI=" << plExponent << ", RMS=" << rms << " (" << diffRMS << ", Q=" << qrms << ", U=" << urms << ", V=" << vrms << "), V%=" << round(10000.0*stokesVFrac)/100.0 << "%" << "(" << sed.AverageFlux(Polarization::StokesV) << "/" << sed.AverageFlux(Polarization::StokesI) << "), max E=" << maxError << " (+" << maxAbsError << ", " << minAbsError << " Jy) @ " << (maxErrorFrequency*1e-6) << ", RMpeak=" << rmPeakValue << " @ " << rmPeakPos << " / " << rmZeroValue << " @0, SB=" << subbandCorrelation << ", 2nd=" << pl2ndOrder;
			return str.str();
		}

		std::string LatexEntry(const Model& model) const
		{
			std::ostringstream str;
			const ModelSource& ms = model.Source(index);
			const MeasuredSED sed = ms.GetIntegratedMSED();
			str
				<< index << " & "
				<< RaDecCoord::RAToString(ms.MeanRA()) << " & "
				<< RaDecCoord::DecToString(ms.MeanDec()) << " & "
				<< round(modalFlux*1000.0) << " mJy & "
				<< round(diffRMS*1000.0) << " mJy & "
				<< round(plExponent*100.0)/100.0;
			return str.str();
		}

	};
	
	std::vector<SourceInfo> _sources;
	
	mutable std::mt19937 _rng;
	
	void make_RA_vs_SI_plot();
	void make_Dec_vs_SI_plot();
	void make_SI_vs_brightness_plot();
	
	void write_csv_file();
	void write_cath_csv_file();
	void outputGeneralSIStats();
	void outputSIBrightnessStatistics(double binStart, double binEnd);
	
	static double sourceRMS(const MeasuredSED& sed, PolarizationEnum pol = Polarization::StokesI);
	
	double sourceRMS(const MeasuredSED& sed, long double factor, long double exponent);
	
	double sourceRMS(const MeasuredSED& sed, long double a, long double b, long double c);
	
	double sourceRMS(const MeasuredSED& sed, const ao::uvector<double>& terms);
	
	static double diffRMS(const MeasuredSED& sed);
	
	double maxStep(const MeasuredSED& sed);
	
	void measureSourceCountNoise();
	
	void getAveragedSED(ao::uvector<double>& averagedValues, const std::set<size_t>& sourceIndices, bool subtractModel, bool useSmoothSpectra);
	void getAveragedSED(MeasuredSED& averagedSED, const std::set<size_t>& sourceIndices, bool subtractModel, bool useSmoothSpectra);
	double getSourceSetWeightSum(const std::set<size_t>& sourceIndices);
	
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

	double getSubbandShapeCorrelation(const MeasuredSED& sed)
	{
		std::vector<double> passband;
		size_t channelsInPassband = 32;
		SubbandPassband::GetSubbandPassband(passband, channelsInPassband);
		makeMeanZero(passband);
		
		double bandStart = sed.begin()->second.FrequencyHz();
		std::vector<double> measurements;
		double correlationSum = 0.0;
		size_t correlationCount = 1;
		for(MeasuredSED::const_iterator m = sed.begin(); m!=sed.end(); ++m)
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
	
	void outputList(const std::vector<SourceInfo>& sortedList, const Model& model, const std::string& name, bool reverse, size_t count);

	void outputSIStats(const std::vector<SourceInfo>& sortedList, size_t n);

	void getAveragedSEDData(ao::uvector<double>& averageResidualFlux, ao::uvector<double>& fluxSigma, size_t nTermsFitted, bool useSmoothSpectra);
	void getAveragedPSData(ao::uvector<double>& power, ao::uvector<double>& powerSigma, size_t nTermsFitted, const class SpectrumFT& sft, bool useSmoothSpectra);
	void runNoiseSimulation(ao::uvector<double>& simPower, size_t nTermsFitted, const MeasuredSED& templateSED, const class SpectrumFT& sft, bool makePsf, bool simulateWithMissingData) const;
	
	void removeEdgeChannels();
	void removeEndChannels();
	void removeBadChannels();
	
	std::string psYAxisDesc()
	{
		if(_psInTemperature)
			return "Power (mK\u00B2[Mpc/h])";
		else
			return "Power (Jy\u00B2)";
	}
};

#endif

