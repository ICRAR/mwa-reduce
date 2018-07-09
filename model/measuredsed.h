#ifndef MEASURED_SED_H
#define MEASURED_SED_H

#include <map>

#include "measurement.h"
#include "spectralenergydistribution.h"

class MeasuredSED : public SpectralEnergyDistribution
{
private:
	typedef std::map<long double, Measurement> FluxMap;
	
public:
	typedef FluxMap::iterator iterator;
	typedef FluxMap::const_iterator const_iterator;
	typedef FluxMap::reverse_iterator reverse_iterator;
	typedef FluxMap::const_reverse_iterator const_reverse_iterator;
	typedef FluxMap::value_type value_type;
		
	MeasuredSED()
	{ }
	
	MeasuredSED(long double fluxDensityJy, long double frequencyHz)
	{
		AddMeasurement(fluxDensityJy, frequencyHz);
	}
	
	template<typename T>
	MeasuredSED(const T* fluxDensityJy, long double frequencyHz)
	{
		AddMeasurement(fluxDensityJy, frequencyHz);
	}
	
	MeasuredSED(long double fluxDensityJy, long double frequencyHz, long double spectralIndex, PolarizationEnum polarization = Polarization::StokesI)
	{
		AddMeasurement(fluxDensityJy, frequencyHz, spectralIndex, polarization);
	}
	
	MeasuredSED(long double fluxDensityAJy, long double frequencyAHz, long double fluxDensityBJy, long double frequencyBHz)
	{
		AddMeasurement(fluxDensityAJy, frequencyAHz);
		AddMeasurement(fluxDensityBJy, frequencyBHz);
	}
	
	MeasuredSED(const MeasuredSED& source)
	: _measurements(source._measurements)
	{
	}
	
	virtual MeasuredSED* Clone() const
	{
		return new MeasuredSED(*this);
	}
	
	void operator=(const MeasuredSED& source)
	{
		_measurements = source._measurements;
	}
	
	void operator+=(const SpectralEnergyDistribution& rhs)
	{
		for(iterator i=begin(); i!=end(); ++i)
		{
			double freq = i->first;
			Measurement &m = i->second;
			for(size_t p=0; p!=4; ++p)
			{
				m.SetFluxDensityFromIndex(p, m.FluxDensityFromIndex(p) + rhs.FluxAtFrequencyFromIndex(freq, p));
			}
		}
	}
	
	void operator*=(double factor)
	{
		for(iterator i=begin(); i!=end(); ++i)
		{
			Measurement &m = i->second;
			for(size_t p=0; p!=4; ++p)
			{
				m.SetFluxDensityFromIndex(p, m.FluxDensityFromIndex(p) * factor);
			}
		}
	}
	
	void CombineMeasurements(const MeasuredSED& other)
	{
		for(const_iterator i=other.begin(); i!=other.end(); ++i)
		{
			double freq = i->first;
			if(_measurements.find(freq) != end())
				throw std::runtime_error("Combining measurements for new frequencies, but frequencies overlap");
			const Measurement& m = i->second;
			AddMeasurement(m);
		}
	}
	
	void CombineMeasurementsWithAveraging(const MeasuredSED& other, double weight = 0.5)
	{
		for(const_iterator i=other.begin(); i!=other.end(); ++i)
		{
			double freq = i->first;
			FluxMap::iterator pos = _measurements.find(freq);
			if(pos == end())
			{
				const Measurement& m = i->second;
				AddMeasurement(m);
			}
			else {
				Measurement& m = pos->second;
				m.AverageWidth(i->second, weight);
			}
		}
	}
	
	void AddMeasurement(const Measurement& measurement)
	{
		_measurements.insert(std::pair<long double, Measurement>(measurement.FrequencyHz(), measurement));
	}
	
	void AddMeasurement(long double fluxDensityJy, long double frequencyHz)
	{
		Measurement measurement;
		measurement.SetFluxDensityFromIndex(0, fluxDensityJy);
		measurement.SetFluxDensityFromIndex(1, 0.0);
		measurement.SetFluxDensityFromIndex(2, 0.0);
		measurement.SetFluxDensityFromIndex(3, 0.0);
		measurement.SetFrequencyHz(frequencyHz);
		_measurements.insert(std::pair<long double, Measurement>(frequencyHz, measurement));
	}
	
	template<typename T>
	void AddMeasurement(const T* fluxDensityJyPerPol, long double frequencyHz)
	{
		Measurement measurement;
		for(size_t p=0; p!=4; ++p)
			measurement.SetFluxDensityFromIndex(p, fluxDensityJyPerPol[p]);
		measurement.SetFrequencyHz(frequencyHz);
		_measurements.insert(std::pair<long double, Measurement>(frequencyHz, measurement));
	}
	
	void AddMeasurement(long double fluxDensityJy, long double frequencyHz, long double spectralIndex, PolarizationEnum polarization = Polarization::StokesI)
	{
#ifdef EXTRA_ASSERTIONS
		if(!Polarization::IsStokes(polarization))
			throw std::runtime_error("Cannot store specified polarization in model");
#endif
		Measurement measurementA, measurementB;
		measurementA.SetZeroExceptSinglePol(polarization, fluxDensityJy);
		measurementA.SetFrequencyHz(frequencyHz);
		_measurements.insert(std::pair<long double, Measurement>(frequencyHz, measurementA));
		if(spectralIndex != 0.0)
		{
			long double fluxB = /* Calculate the flux density for 1 Hz frequency */
				fluxDensityJy * std::pow(1.0 / frequencyHz, spectralIndex);
			long double refFreqB = frequencyHz + 15000000.0;
			if(refFreqB == frequencyHz) {
				refFreqB *= 2.0;
			}
			fluxB = fluxB * std::pow(refFreqB, spectralIndex);
			measurementB.SetZeroExceptSinglePol(polarization, fluxB);
			measurementB.SetFrequencyHz(refFreqB);
			_measurements.insert(std::pair<long double, Measurement>(refFreqB, measurementB));
		}
	}
	
	std::string ToString() const
	{
		std::ostringstream s;
		s.precision(15);
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			i->second.ToStream(s);
		}
		return s.str();
	}
	
	long double FluxAtFrequencyFromIndex(long double frequencyHz, size_t pIndex) const
	{
		if(_measurements.size() <= 1)
		{
			if(_measurements.empty()) return 0.0;
			else return _measurements.begin()->second.FluxDensityFromIndex(pIndex);
		}
		
		// 'right' will be first item which frequency >= frequencyHz
		FluxMap::const_iterator right = _measurements.lower_bound(frequencyHz);
		if(right != _measurements.end() && right->first == frequencyHz)
			return right->second.FluxDensityFromIndex(pIndex);
		
		FluxMap::const_iterator left;
		
		// If the requested frequency is outside the range, we extrapolate the SI of the full range
		if(right == _measurements.begin() || right == _measurements.end())
		{
			left = _measurements.begin();
			right = _measurements.end();
			--right;
		} else {
			// Requested frequency is within range (no extrapolation required)
			left = right;
			--left;
		}
		
		long double
			freqA = left->first,
			fluxA = left->second.FluxDensityFromIndex(pIndex),
			freqB = right->first,
			fluxB = right->second.FluxDensityFromIndex(pIndex);
		
		return SpectralEnergyDistribution::FluxAtFrequency(fluxA, freqA, fluxB, freqB, frequencyHz);
	}
	
	long double IntegratedFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		if(startFrequency == endFrequency)
			return FluxAtFrequency(startFrequency, polarization);
		
		FluxMap::const_iterator iter = _measurements.lower_bound(startFrequency);
		
		/** Handle special cases */
		if(_measurements.size() <= 2)
		{
			if(_measurements.empty())
				return 0.0;
			else if(_measurements.size()==1)
				return _measurements.begin()->second.FluxDensity(polarization);
			else { // _measurements.size()==2
				long double
					freqA = _measurements.begin()->first,
					fluxA = _measurements.begin()->second.FluxDensity(polarization),
					freqB = _measurements.rbegin()->first,
					fluxB = _measurements.rbegin()->second.FluxDensity(polarization);
				return SpectralEnergyDistribution::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
		}
		if(iter == _measurements.end()) { // all keys are lower, so take entire range
			long double
				freqA = _measurements.begin()->first,
				fluxA = _measurements.begin()->second.FluxDensity(polarization),
				freqB = _measurements.rbegin()->first,
				fluxB = _measurements.rbegin()->second.FluxDensity(polarization);
			return SpectralEnergyDistribution::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
		}
		
		if(iter != _measurements.begin()) --iter;
		
		if(iter->first >= endFrequency) {
			// all keys are outside range, higher than range
			long double
				freqA = _measurements.begin()->first,
				fluxA = _measurements.begin()->second.FluxDensity(polarization),
				freqB = _measurements.rbegin()->first,
				fluxB = _measurements.rbegin()->second.FluxDensity(polarization);
			return SpectralEnergyDistribution::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
		}
		
		long double integratedSum = 0.0;
		long double leftFrequency = startFrequency;
		if(leftFrequency < iter->first)
		{
			// requested frequency is below first item; extrapolate
			long double
				freqA = _measurements.begin()->first,
				fluxA = _measurements.begin()->second.FluxDensity(polarization),
				freqB = _measurements.rbegin()->first,
				fluxB = _measurements.rbegin()->second.FluxDensity(polarization);
			long double sumTerm = SpectralEnergyDistribution::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, iter->first);
			integratedSum += sumTerm * (iter->first - startFrequency);
			leftFrequency = iter->first;
		}
			
		while(iter != _measurements.end() && iter->first < endFrequency)
		{
			FluxMap::const_iterator left = iter;
			FluxMap::const_iterator right = iter;
			++right;
			
			long double rightFrequency;
			
			// If this is past the sampled frequencies, extrapolate full range
			if(right == _measurements.end()) {
				left = _measurements.begin();
				right = _measurements.end();
				--right;
				rightFrequency = endFrequency;
			} else {
				rightFrequency = right->first;
				if(rightFrequency > endFrequency)
					rightFrequency = endFrequency;
			}
			
			long double
				freqA = left->first,
				fluxA = left->second.FluxDensity(polarization),
				freqB = right->first,
				fluxB = right->second.FluxDensity(polarization);
			
			if(leftFrequency < rightFrequency)
			{
				long double sumTerm = SpectralEnergyDistribution::IntegratedFlux(fluxA, freqA, fluxB, freqB, leftFrequency, rightFrequency);
				if(!std::isfinite(sumTerm))
				{
					std::cerr << "Warning: integrating flux between " << leftFrequency << " and " << rightFrequency << " with fluxes " << fluxA << '@' << freqA << ',' << fluxB << '@' << freqB << " gave non-finite result\n";
				}
				
				integratedSum += sumTerm * (rightFrequency - leftFrequency);
			}
			leftFrequency = rightFrequency;
			++iter;
		}
		return integratedSum / (endFrequency - startFrequency);
	}
	
	long double AverageFlux(PolarizationEnum polarization) const
	{
		long double sum = 0.0;
		size_t count = 0;
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			const Measurement &m = i->second;
			long double flux = m.FluxDensity(polarization);
			if(std::isfinite(flux))
			{
				++count;
				sum += flux;
			}
		}
		return sum / (long double) count;
	}
	
	long double AverageFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const
	{
		if(startFrequency == endFrequency)
			return FluxAtFrequency(startFrequency, polarization);
		
		/** Handle special cases */
		if(_measurements.empty())
			return 0.0;
		
		FluxMap::const_iterator iter = _measurements.lower_bound(startFrequency);
		if(iter == _measurements.end()) // all keys are lower
			return std::numeric_limits<long double>::quiet_NaN();
		
		size_t count = 0;
		long double fluxSum = 0.0;
		
		while(iter->first < endFrequency && iter != _measurements.end())
		{
			long double flux = iter->second.FluxDensity(polarization);
			if(std::isfinite(flux))
			{
				++count;
				fluxSum += flux;
			}
			++iter;
		}
		
		if(count == 0)
			return std::numeric_limits<long double>::quiet_NaN();
		else
			return fluxSum / count;
	}
	
	void FitPowerlaw(long double& factor, long double& exponent, PolarizationEnum polarization) const
	{
		long double sumxy = 0.0, sumx = 0.0, sumy = 0.0, sumxx = 0.0;
		size_t n = 0;
		bool requireNonLinear = false;
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			const Measurement &m = i->second;
			long double flux = m.FluxDensity(polarization);
			if(flux <= 0)
			{
				requireNonLinear = true;
				break;
			}
			if(m.FrequencyHz() > 0 && flux > 0 && std::isfinite(flux))
			{
				long double
					logx = std::log(m.FrequencyHz()),
					logy = std::log(flux);
				sumxy += logx * logy;
				sumx += logx;
				sumy += logy;
				sumxx += logx * logx;
				++n;
			}
		}
		if(requireNonLinear)
		{
			NonLinearPowerLawFitter fitter;
			n = 0;
			for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
			{
				const Measurement &m = i->second;
				long double flux = m.FluxDensity(polarization);
				if(std::isfinite(m.FrequencyHz()) && std::isfinite(flux)) {
					fitter.AddDataPoint(m.FrequencyHz(), flux);
					++n;
				}
			}
			double eTemp = 0.0, fTemp = 1.0;
			fitter.Fit(eTemp, fTemp);
			//if(n == 0)
			//	std::cout << "No valid data in power law fit\n";
			//else
			//	std::cout << "Non-linear fit yielded: " << fTemp << " * x^" << eTemp << "\n";
			exponent = eTemp;
			factor = fTemp;
		}
		else {
			if(n == 0)
			{
				exponent = std::numeric_limits<double>::quiet_NaN();
				factor = std::numeric_limits<double>::quiet_NaN();
			}
			else {
				exponent = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
				factor = std::exp((sumy - exponent * sumx) / n);
			}
		}
	}
	
	void FitPowerlaw2ndOrder(long double& a, long double& b, long double& c, PolarizationEnum polarization) const
	{
		NonLinearPowerLawFitter fitter;
		size_t n = 0;
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			const Measurement &m = i->second;
			long double flux = m.FluxDensity(polarization);
			if(std::isfinite(m.FrequencyHz()) && std::isfinite(flux)) {
				fitter.AddDataPoint(m.FrequencyHz(), flux);
				++n;
			}
		}
		double aTemp = 0.0, bTemp = 1.0, cTemp = 0.0;
		fitter.Fit(aTemp, bTemp, cTemp);
		a = aTemp;
		b = bTemp;
		c = cTemp;
	}
	
	
	void FitLogPolynomial(ao::uvector<double>& terms, size_t nTerms, PolarizationEnum polarization, double referenceFrequencyHz = 1.0) const
	{
		NonLinearPowerLawFitter fitter;
		size_t n = 0;
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			const Measurement &m = i->second;
			long double flux = m.FluxDensity(polarization);
			if(std::isfinite(m.FrequencyHz()) && std::isfinite(flux)) {
				fitter.AddDataPoint(m.FrequencyHz() / referenceFrequencyHz, flux);
				++n;
			}
		}
		terms.assign(nTerms, 0.0);
		fitter.Fit(terms, nTerms);
	}
	
	virtual long double ReferenceFrequencyHz() const
	{
		return (_measurements.begin()->second.FrequencyHz() +
			_measurements.rbegin()->second.FrequencyHz()) * 0.5;
	}
	
	long double FluxAtLowestFrequency() const
	{
		const Measurement &m = _measurements.begin()->second;
		return m.FluxDensityFromIndex(0);
	}
	
	bool HasValidMeasurement() const
	{
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			if(std::isfinite(i->second.FluxDensity(Polarization::StokesI)) ||
				std::isfinite(i->second.FluxDensity(Polarization::StokesQ)) ||
				std::isfinite(i->second.FluxDensity(Polarization::StokesU)) ||
				std::isfinite(i->second.FluxDensity(Polarization::StokesV)))
				return true;
		}
		return false;
	}
	
	void RemoveInvalidMeasurements()
	{
		FluxMap::iterator i=_measurements.begin();
		while(i!=_measurements.end())
		{
			if(!std::isfinite(i->second.FluxDensity(Polarization::StokesI)) ||
				!std::isfinite(i->second.FluxDensity(Polarization::StokesQ)) ||
				!std::isfinite(i->second.FluxDensity(Polarization::StokesU)) ||
				!std::isfinite(i->second.FluxDensity(Polarization::StokesV)))
			{
				_measurements.erase(i);
				i = _measurements.begin();
			}
			else {
				++i;
			}
		}
	}
	
	bool operator<(const SpectralEnergyDistribution& other) const
	{
		const MeasuredSED& msed = dynamic_cast<const MeasuredSED&>(other);
		double thisFrequency = _measurements.begin()->first;
		double otherFrequency = msed._measurements.begin()->first;
		double minFreq = std::min(thisFrequency, otherFrequency);
		return
			FluxAtFrequencyFromIndex(minFreq, 0) < msed.FluxAtFrequencyFromIndex(minFreq, 0);
	}
	
	size_t MeasurementCount() const { return _measurements.size(); }
	long double LowestFrequency() const { return _measurements.begin()->first; }
	long double HighestFrequency() const { return _measurements.rbegin()->first; }
	long double CentreFrequency() const { return 0.5 * (LowestFrequency() + HighestFrequency()); }
	
	void GetMeasurements(std::vector<Measurement> &measurements) const
	{
		for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
		{
			measurements.push_back(i->second);
		}
	}
	
	//iterator begin() { return boost::adaptors::values(_measurements).begin(); }
	iterator begin() { return _measurements.begin(); }
	const_iterator begin() const { return _measurements.begin(); }
	iterator end() { return _measurements.end(); }
	const_iterator end() const { return _measurements.end(); }
	
	reverse_iterator rbegin() { return _measurements.rbegin(); }
	const_reverse_iterator rbegin() const { return _measurements.rbegin(); }
	reverse_iterator rend() { return _measurements.rend(); }
	const_reverse_iterator rend() const { return _measurements.rend(); }
private:
	FluxMap _measurements;
};

#endif
