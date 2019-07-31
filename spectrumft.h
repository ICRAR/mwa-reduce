#ifndef SPECTRUM_FT_H
#define SPECTRUM_FT_H

#include "uvector.h"
#include "model/measuredsed.h"
#include "universe.h"
#include "units/angle.h"

#include <fftw3.h>

#include <cstring>
#include <complex>

class SpectrumFT
{
public:
	explicit SpectrumFT(size_t channelCount);
	~SpectrumFT();
	
	void GetFTPower(ao::uvector<double>& destination, const MeasuredSED& sed, size_t nTermsFitted, bool toTemperature, bool useLomb) const;
	
	double GetAveragePowerInKRange(double kLow, double kHigh, const MeasuredSED& sed, size_t nTermsFitted, bool toTemperature, bool useLomb);
	
	static double GetMaxKParallell(const MeasuredSED& sed, bool verbose=true);
	
	static double GetMaxDelayInMuSec(const MeasuredSED& sed);
	
	static double FrequencyToDistance(double frequency);
	
	static void Normalize(ao::uvector<double>& ps)
	{
		double f = 0.0;
		for(ao::uvector<double>::iterator i=ps.begin(); i!=ps.end(); ++i)
			f += (*i) > 0.0 ? (*i) : 0.0;
		f = double(ps.size()) / f;
		for(ao::uvector<double>::iterator i=ps.begin(); i!=ps.end(); ++i)
			(*i) *= f;
	}
	static double LowerKParallellEoRWindow()
	{
		return 6.5e-2;
	}
	static double UpperKParallellEoRWindow()
	{
		return 3.7;
	}
	
	static void ConvertPSToTemperature(ao::uvector<double>& ps, const double* frequencies, const double* pBeamOmegas)
	{
		// Good references:
		// - Furlanetto et al 2006
		// - http://arxiv.org/abs/1103.2135 Parsons et al 2012 ApJ
		double bandwidth = frequencies[ps.size()-1] - frequencies[0];
		double centralOmega = pBeamOmegas[ps.size()/2];
		double centralFreq = frequencies[0] + bandwidth*0.5;
		const double c = Universe::C(), b = Universe::Boltzmann();
		double distance = FrequencyToDistance(centralFreq);
		double nuToR = distance / centralFreq; // parallel
		double lToR = distance;                // transverse
		std::cout << "distance=" << distance << ", nuToR=" << nuToR << ", lToR=" << lToR << '\n';
		
		double f = frequencies[ps.size()/2];
		double term = (2.0 * b * f * f) / (c*c);
		double factor = /*flux to watts*/ 1e-23 * /*K^2 to mK^2*/ 1e-6 * (lToR * lToR * nuToR) / (term * term * bandwidth * centralOmega);
			
		for(size_t channelIndex=0; channelIndex!=ps.size(); ++channelIndex)
		{
			if(channelIndex == ps.size()/2)
				std::cout << "Channel " << channelIndex << " ps to temp factor = " << factor << '\n';
			ps[channelIndex] *= factor;
		}
	}
	
	void ConvertPSToTemperature(ao::uvector<double>& ps)
	{
		ConvertPSToTemperature(ps, _frequencies.data(), _pBeamOmegas.data());
	}
	
	static void ConvertFluxToTemperature(ao::uvector<double>& data, const double* frequencies, const double* sBeamOmegas)
	{
		const double c = Universe::C(), b = Universe::Boltzmann();
		for(size_t channelIndex=0; channelIndex!=data.size(); ++channelIndex)
		{
			double f = frequencies[channelIndex];
			//double lambda = c/f;
			// conversion to mK
			double jToT = /*2.0*M_PI **/ 1e-23 * c*c / (2.0 * b * f * f * sBeamOmegas[channelIndex]);
			data[channelIndex] *= jToT;
		}
	}
	
	void ConvertFluxToTemperature(ao::uvector<double>& data)
	{
		ConvertFluxToTemperature(data, _frequencies.data(), _sBeamOmegas.data());
	}
	
	void SetMetaData(const MeasuredSED& sed, double maximumBaselineInM, size_t sourceCount)
	{
		ao::uvector<double>
			frequencies(sed.MeasurementCount());
		MeasuredSED::const_iterator sedIter = sed.begin();
		for(size_t i=0; i!=sed.MeasurementCount(); ++i)
		{
			frequencies[i] = sedIter->first;
			++sedIter;
		}
		SetMetaData(frequencies, maximumBaselineInM, sourceCount);
	}
	
	void SetMetaData(const ao::uvector<double>& frequencies, double maximumBaselineInM, size_t sourceCount)
	{
		_frequencies = frequencies;
		_sBeamOmegas.resize(_frequencies.size());
		_pBeamOmegas.resize(_frequencies.size());
		const double c = Universe::C();
		for(size_t i=0; i!=frequencies.size(); ++i)
		{
			double f = frequencies[i];
			
			double maxBaselineInLambda = maximumBaselineInM*f/c;
			double sBeam = 1.0 / maxBaselineInLambda;
			_sBeamOmegas[i] = M_PI*sBeam*sBeam/(4.0*M_LN2);
			
			//double elementSizeInLambda = elementSizeInM*f/c;
			//double pBeam = 1.0 / elementSizeInLambda;
			//_pBeamOmegas[i] = M_PI*pBeam*pBeam/(4.0*M_LN2);
			_pBeamOmegas[i] = _sBeamOmegas[i] * sourceCount;
			
			if(i==frequencies.size()/2)
			{
				std::cout << "Synthesized beam at " << f*1e-6 << " MHz: " << Angle::ToNiceString(sBeam) << " FWHM (omega=" << _sBeamOmegas[i] << ")\n";
				double pb = sqrt(_pBeamOmegas[i] / (M_PI/(4.0*M_LN2)));
				std::cout << "Effective primary beam at " << f*1e-6 << " MHz: " << Angle::ToNiceString(pb) << " FWHM (omega=" << _pBeamOmegas[i] << ")\n";
			}
		}
	}
	void SetOutputFit(bool outputFit) {
		_outputFit = outputFit;
	}
	void PreparePSPlot(class GNUPlot& plot, double centralFrequency, bool useDelay);
private:
	void transformSpectrum(ao::uvector<std::complex<double>>& data, size_t sampleCount) const;
	void lombPeriodogram(ao::uvector<double>& output, const ao::uvector<double>& input) const;
	void twoPoint(ao::uvector<double>& output, const ao::uvector<double>& input) const;
	void applyWindow(ao::uvector<double>& data) const;
	static double getKForFreqStep(double freqStepHz, double centralFreq);
	void interpolate(ao::uvector<double>& data, size_t i) const;
	void interpolate(ao::uvector<double>& data) const;
	
	size_t _channelCount;
	fftw_plan _fftPlan;
	bool _outputFit;
	ao::uvector<double> _frequencies;
	ao::uvector<double> _sBeamOmegas, _pBeamOmegas;
	std::complex<double>* _inputData;
	std::complex<double>* _outputData;
};

#endif
