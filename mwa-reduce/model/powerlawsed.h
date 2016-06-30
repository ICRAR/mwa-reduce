#ifndef POWER_LAW_SED_H
#define POWER_LAW_SED_H

#include "spectralenergydistribution.h"

class PowerLawSED : public SpectralEnergyDistribution
{
public:
	PowerLawSED() : _referenceFrequency(0.0)
	{
		for(size_t p=0; p!=4; ++p) _factors[p] = 0.0;
	}
	
	PowerLawSED(double referenceFrequency, double constantFlux) :
		_referenceFrequency(referenceFrequency),
		_terms(1)
	{
		double refBrightness = constantFlux;
		if(refBrightness <= 0.0)
			refBrightness = 1.0;
		_terms[0] = refBrightness;
		_factors[0] = constantFlux / refBrightness;
		for(size_t p=1; p!=4; ++p)
			_factors[p] = 0.0;
	}
	
	virtual PowerLawSED* Clone() const
	{
		return new PowerLawSED(*this);
	}
	
	virtual std::string ToString() const
	{
		std::ostringstream str;
		double term1 = _terms.size()>1 ? _terms[1] : 0.0;
		double f = _terms[0];
		double
			i=f*_factors[0], q=f*_factors[1],
			u=f*_factors[2], v=f*_factors[3];
		str << "    sed {\n      frequency " << _referenceFrequency*1e-6 << " MHz\n      fluxdensity Jy "
			<< i << " " << q << " " << u << " " << v <<
			"\n      spectral-index { ";
		str << term1;
		for(size_t i=2; i<_terms.size(); ++i)
			str << ", " << _terms[i];
		str << " }\n    }\n";
		return str.str();
	}
	
	virtual long double FluxAtFrequencyFromIndex(long double frequencyHz, size_t pIndex) const
	{
		return NonLinearPowerLawFitter::Evaluate(frequencyHz/_referenceFrequency, _terms) * _factors[pIndex];
	}
	
	virtual long double IntegratedFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		size_t pIndex = Polarization::StokesToIndex(polarization);
		long double sum = 0.0;
		for(size_t i=0; i!=101; ++i)
		{
			long double frequency = startFrequency + (endFrequency - startFrequency)*double(i)/100.0;
			sum += FluxAtFrequencyFromIndex(frequency, pIndex);
		}
		return sum / 101.0;
	}
	
	virtual long double AverageFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		return IntegratedFlux(startFrequency, endFrequency, polarization);
	}
	
	virtual bool operator<(const SpectralEnergyDistribution &other) const
	{
		return other.FluxAtFrequencyFromIndex(_referenceFrequency, 0) < FluxAtFrequencyFromIndex(_referenceFrequency, 0);
	}
	
	virtual void operator*=(double factor)
	{
		for(size_t p=0; p!=4; ++p)
			_factors[p] *= factor;
	}
	
	virtual void operator+=(const SpectralEnergyDistribution &other)
	{
		throw std::runtime_error("operator+= not yet implemented for power law sed");
	}
	
	virtual long double ReferenceFrequencyHz() const
	{
		return _referenceFrequency;
	}
	
	void SetData(double referenceFrequency, const double* brightnessVector, const ao::uvector<double>& siTerms)
	{
		_referenceFrequency = referenceFrequency;
		double refBrightness = brightnessVector[0];
		if(refBrightness <= 0.0)
			refBrightness = 1.0;
		_terms.resize(siTerms.size()+1);
		_terms[0] = refBrightness;
		for(size_t i=0; i!=siTerms.size(); ++i)
			_terms[i+1] = siTerms[i];
		for(size_t p=0; p!=4; ++p)
			_factors[p] = brightnessVector[p] / refBrightness;
	}
	
	void GetData(double& referenceFrequency, double* brightnessVector, std::vector<double>& siTerms) const
	{
		referenceFrequency = _referenceFrequency;
		double f = _terms[0];
		for(size_t p=0; p!=4; ++p)
			brightnessVector[p] = f*_factors[p];
		siTerms.resize(_terms.size()-1);
		for(size_t i=0; i!=_terms.size()-1; ++i)
			siTerms[i] = _terms[i+1];
	}
private:
	double _referenceFrequency;
	double _factors[4];
	ao::uvector<double> _terms;
};

#endif
