#ifndef SPECTRAL_ENERGY_DISTRIBUTION_H
#define SPECTRAL_ENERGY_DISTRIBUTION_H

#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

class Measurement
{
	public:
		Measurement() :
			_isApparent(false),
			_frequencyHz(0.0),
			_bandWidthHz(0.0)
		{
		}
		
		Measurement(const Measurement &source) :
			_isApparent(source._isApparent),
			_frequencyHz(source._frequencyHz),
			_bandWidthHz(source._bandWidthHz)
		{
			memcpy(_fluxDensities, source._fluxDensities, sizeof(long double)*4);
			memcpy(_fluxDensityStddevs, source._fluxDensityStddevs, sizeof(long double)*4);
			memcpy(_beamValue, source._beamValue, sizeof(long double)*4);
		}
		
		void operator=(const Measurement &source)
		{ 
			_isApparent = source._isApparent;
			_frequencyHz = source._frequencyHz;
			_bandWidthHz = source._bandWidthHz;
			memcpy(_fluxDensities, source._fluxDensities, sizeof(long double)*4);
			memcpy(_fluxDensityStddevs, source._fluxDensityStddevs, sizeof(long double)*4);
			memcpy(_beamValue, source._beamValue, sizeof(long double)*4);
		}
		
		long double FrequencyHz() const { return _frequencyHz; }
		
		void SetFrequencyHz(long double frequencyHz) { _frequencyHz = frequencyHz; }
		
		long double FluxDensity(size_t polarizationIndex) const { return _fluxDensities[polarizationIndex]; }
		
		void SetFluxDensity(size_t polarizationIndex, long double flux) { _fluxDensities[polarizationIndex] = flux; }
		
		void SetFluxDensityStddev(size_t polarizationIndex, long double stddev) { _fluxDensityStddevs[polarizationIndex] = stddev; }
		
		void SetIsApparent(bool isApparent) { _isApparent = isApparent; }
		
		void SetBandWidthHz(double bandwidthHz) { _bandWidthHz = bandwidthHz; }
		
		void SetBeamValue(size_t polarizationIndex, long double flux) { _beamValue[polarizationIndex] = flux; }
		
		void ToStream(std::ostream &s) const
		{
			s <<
				"    measurement {\n"
				"      type ";
			if(_isApparent)
				s << "apparent";
			else
				s << "absolute";
			s << "\n"
				"      frequency " << (_frequencyHz/1000000.0) << " MHz\n"
				"      fluxdensity Jy " << _fluxDensities[0] << ' ' << _fluxDensities[1] << ' '
				<< _fluxDensities[2] << ' ' << _fluxDensities[3] << '\n';
			if(_bandWidthHz > 0.0)
				s << "      bandwidth " << _bandWidthHz << " Hz\n";
			if(_beamValue[0] > 0.0 || _beamValue[1] > 0.0 ||_beamValue[2] > 0.0 ||_beamValue[3] > 0.0)
				s << "      beam-value " << _beamValue[0] << ' ' << _beamValue[1] << ' '
				<< _beamValue[2] << ' ' << _beamValue[3] << '\n';
			s << "    }\n";
		}
	private:
		
		bool _isApparent;
		double _frequencyHz;
		double _bandWidthHz;
		long double _fluxDensities[4];
		long double _fluxDensityStddevs[4];
		long double _beamValue[4];
};

class SpectralEnergyDistribution
{
	private:
		typedef std::map<long double, Measurement> FluxMap;
	public:
		SpectralEnergyDistribution()
		{
		}
		
		SpectralEnergyDistribution(long double fluxDensityJy, long double frequencyHz)
		{
			AddMeasurement(fluxDensityJy, frequencyHz);
		}
		
		SpectralEnergyDistribution(long double fluxDensityJy, long double frequencyHz, long double spectralIndex)
		{
			AddMeasurement(fluxDensityJy, frequencyHz, spectralIndex);
		}
		
		SpectralEnergyDistribution(long double fluxDensityAJy, long double frequencyAHz, long double fluxDensityBJy, long double frequencyBHz)
		{
			AddMeasurement(fluxDensityAJy, frequencyAHz);
			AddMeasurement(fluxDensityBJy, frequencyBHz);
		}
		
		SpectralEnergyDistribution(const SpectralEnergyDistribution &source)
		: _measurements(source._measurements)
		{
		}
		
		void operator=(const SpectralEnergyDistribution &source)
		{
			_measurements = source._measurements;
		}
		
		void AddMeasurement(const Measurement &measurement)
		{
			_measurements.insert(std::pair<long double, Measurement>(measurement.FrequencyHz(), measurement));
		}
		
		void AddMeasurement(long double fluxDensityJy, long double frequencyHz)
		{
			Measurement measurement;
			measurement.SetFluxDensity(0, fluxDensityJy);
			measurement.SetFluxDensity(1, 0.0);
			measurement.SetFluxDensity(2, 0.0);
			measurement.SetFluxDensity(3, fluxDensityJy);
			measurement.SetFrequencyHz(frequencyHz);
			_measurements.insert(std::pair<long double, Measurement>(frequencyHz, measurement));
		}
		
		void AddMeasurement(long double fluxDensityJy, long double frequencyHz, long double spectralIndex)
		{
			Measurement measurementA, measurementB;
			measurementA.SetFluxDensity(0, fluxDensityJy);
			measurementA.SetFluxDensity(1, 0.0);
			measurementA.SetFluxDensity(2, 0.0);
			measurementA.SetFluxDensity(3, fluxDensityJy);
			measurementA.SetFrequencyHz(frequencyHz);
			_measurements.insert(std::pair<long double, Measurement>(frequencyHz, measurementA));
			if(spectralIndex != 0.0)
			{
				long double fluxB = /* Calculate the flux density for 1 Hz frequency */
					fluxDensityJy * std::pow(1.0 / frequencyHz, spectralIndex);
				long double refFreqB = 1.0;
				if(refFreqB == frequencyHz) {
					refFreqB = 2.0;
					fluxB = fluxB * std::pow(frequencyHz, spectralIndex);
				}
				measurementB.SetFluxDensity(0, fluxB);
				measurementB.SetFluxDensity(1, 0.0);
				measurementB.SetFluxDensity(2, 0.0);
				measurementB.SetFluxDensity(3, fluxB);
				measurementB.SetFrequencyHz(refFreqB);
				_measurements.insert(std::pair<long double, Measurement>(refFreqB, measurementB));
			}
		}
		
		void SetConstantBeam(double xx, double xy, double yx, double yy)
		{
			for(std::map<long double, Measurement>::iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
			{
				i->second.SetBeamValue(0, xx);
				i->second.SetBeamValue(1, xy);
				i->second.SetBeamValue(2, yx);
				i->second.SetBeamValue(3, yy);
			}
		}
		
		std::string ToString() const
		{
			std::ostringstream s;
			for(FluxMap::const_iterator i=_measurements.begin(); i!=_measurements.end(); ++i)
			{
				i->second.ToStream(s);
			}
			return s.str();
		}
		
		long double FluxAtFrequency(long double frequencyHz, size_t polarizationIndex) const
		{
			if(_measurements.size() <= 1)
			{
				if(_measurements.empty()) return 0.0;
				else return _measurements.begin()->second.FluxDensity(polarizationIndex);
			}
			
			// 'right' will be first item which frequency >= frequencyHz
			FluxMap::const_iterator right = _measurements.lower_bound(frequencyHz);
			if(right->first == frequencyHz)
				return right->second.FluxDensity(polarizationIndex);
			
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
				fluxA = left->second.FluxDensity(polarizationIndex),
				freqB = right->first,
				fluxB = right->second.FluxDensity(polarizationIndex);
			
			return FluxAtFrequency(fluxA, freqA, fluxB, freqB, frequencyHz);
		}
		
		static long double FluxAtFrequency(long double fluxDensityAJy, long double referenceFrequencyAHz,
																			 long double fluxDensityBJy, long double referenceFrequencyBHz,
																			long double requestedFrequency)
		{
			// if either fluxes are zero, or one of them is negative and the other not,
			// perform linear interpolation instead of power law interpolation
			bool signA = fluxDensityAJy < 0.0, signB = fluxDensityBJy < 0.0;
			if(fluxDensityAJy==0.0 || fluxDensityBJy==0.0 || (signA && !signB) || (signB && !signA))
			{
				long double slope =
					(fluxDensityBJy - fluxDensityAJy) /
					(referenceFrequencyBHz - referenceFrequencyAHz);
				return fluxDensityAJy + slope * (requestedFrequency - referenceFrequencyAHz);
			} else {
				long double si =
					log(fabs(fluxDensityBJy/fluxDensityAJy)) /
					log(referenceFrequencyBHz/referenceFrequencyAHz);
				return fluxDensityAJy * std::pow(requestedFrequency/referenceFrequencyAHz, si);
			}
		}
		
		long double FluxAtChannel(size_t channelIndex, size_t channelCount, long double startFreq, long double endFreq, size_t polarizationIndex) const
		{
			long double freq = startFreq + (long double) channelIndex * (endFreq - startFreq) / (long double) (channelCount-1);
			return FluxAtFrequency(freq, polarizationIndex);
		}
		
		long double IntegratedFlux(long double startFrequency, long double endFrequency, size_t polarizationIndex) const
		{
			if(startFrequency == endFrequency)
				return FluxAtFrequency(startFrequency, polarizationIndex);
			
			FluxMap::const_iterator iter = _measurements.lower_bound(startFrequency);
			
			/** Handle special cases */
			if(_measurements.size() <= 2)
			{
				if(_measurements.empty())
					return 0.0;
				else if(_measurements.size()==1)
					return _measurements.begin()->second.FluxDensity(polarizationIndex);
				else if(_measurements.size()==2) {
					long double
						freqA = _measurements.begin()->first,
						fluxA = _measurements.begin()->second.FluxDensity(polarizationIndex),
						freqB = _measurements.rbegin()->first,
						fluxB = _measurements.rbegin()->second.FluxDensity(polarizationIndex);
					
					return IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
				}
			}
			if(iter == _measurements.end()) { // all keys are lower, so take entire range
				long double
					freqA = _measurements.begin()->first,
					fluxA = _measurements.begin()->second.FluxDensity(polarizationIndex),
					freqB = _measurements.rbegin()->first,
					fluxB = _measurements.rbegin()->second.FluxDensity(polarizationIndex);
				return IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
			
			if(iter != _measurements.begin()) --iter;
			
			if(iter->first >= endFrequency) {
				// all keys are outside range, higher than range
				long double
					freqA = _measurements.begin()->first,
					fluxA = _measurements.begin()->second.FluxDensity(polarizationIndex),
					freqB = _measurements.rbegin()->first,
					fluxB = _measurements.rbegin()->second.FluxDensity(polarizationIndex);
				return IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
			
			long double integratedSum = 0.0;
			long double leftFrequency = startFrequency;
			if(leftFrequency < iter->first)
			{
				// requested frequency is below first item; extrapolate
				long double
					freqA = _measurements.begin()->first,
					fluxA = _measurements.begin()->second.FluxDensity(polarizationIndex),
					freqB = _measurements.rbegin()->first,
					fluxB = _measurements.rbegin()->second.FluxDensity(polarizationIndex);
				long double sumTerm = IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, iter->first);
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
					fluxA = left->second.FluxDensity(polarizationIndex),
					freqB = right->first,
					fluxB = right->second.FluxDensity(polarizationIndex);
				
				if(leftFrequency < rightFrequency)
				{
					long double sumTerm = IntegratedFlux(fluxA, freqA, fluxB, freqB, leftFrequency, rightFrequency);
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
		
		static long double IntegratedFlux(long double fluxDensityAJy, long double referenceFrequencyAHz, long double fluxDensityBJy, long double referenceFrequencyBHz, long double startFrequency, long double endFrequency)
		{
			// if either fluxes are zero, or one of them is negative and the other not,
			// perform linear interpolation instead of power law interpolation
			bool signA = fluxDensityAJy < 0.0, signB = fluxDensityBJy < 0.0;
			if(fluxDensityAJy==0.0 || fluxDensityBJy==0.0 || (signA && !signB) || (signB && !signA))
			{
				long double slope =
					(fluxDensityBJy - fluxDensityAJy) /
					(referenceFrequencyBHz - referenceFrequencyAHz);
				return slope * 0.5 * (
					endFrequency*endFrequency - startFrequency*startFrequency) /
					(endFrequency - startFrequency)
					+ (fluxDensityAJy - slope * referenceFrequencyAHz);
			} else {
				long double si =
					log(fluxDensityBJy/fluxDensityAJy) /
					log(referenceFrequencyBHz/referenceFrequencyAHz);
				return fluxDensityAJy * (
					std::pow(endFrequency/referenceFrequencyAHz, si) * endFrequency -
					std::pow(startFrequency/referenceFrequencyAHz, si) * startFrequency) /
					( (si+1.0) * (endFrequency - startFrequency) );
			}
		}
		
	private:
		FluxMap _measurements;
};

#endif
