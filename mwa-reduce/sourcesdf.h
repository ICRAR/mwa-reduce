#ifndef SOURCE_STRENGTH_H
#define SOURCE_STRENGTH_H

#include <cmath>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

template<typename NumericType=long double>
class SourceSDF
{
	public:
		virtual ~SourceSDF() { }
		
		virtual NumericType FluxAtFrequency(NumericType frequencyHz) const = 0;
			
		virtual NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const = 0;
		
		NumericType FluxAtFrequency(size_t channelIndex, size_t channelCount, NumericType startFreq, NumericType endFreq) const
		{
			NumericType freq = startFreq + NumericType(channelIndex) * (endFreq - startFreq) / NumericType(channelCount-1);
			return FluxAtFrequency(freq);
		}
		
		virtual SourceSDF<NumericType> *Copy() const = 0;
		
		virtual std::string ToString() const = 0;
		
		static SourceSDF<NumericType> *ParseLine(std::vector<std::string>::const_iterator firstToken, std::vector<std::string>::const_iterator endToken);
};

template<typename NumericType=long double>
class SourceSDFWithSI : public SourceSDF<NumericType>
{
	public:
		SourceSDFWithSI() :
		_fluxDensityJy(0.0),
		_spectralIndex(0.0),
		_refFreqA(1e8),
		_refFreqB(2e8)
		{
		}
		
		SourceSDFWithSI(NumericType fluxDensityJy, NumericType spectralIndex, NumericType siReferenceFrequencyHz) :
		/* Calculate the flux density for 1 Hz frequency */
			_fluxDensityJy( fluxDensityJy *
				std::pow((NumericType) 1.0 / siReferenceFrequencyHz, spectralIndex) ),
			_spectralIndex( spectralIndex ),
			_refFreqA(siReferenceFrequencyHz),
			_refFreqB(siReferenceFrequencyHz*2.0)
		{
		}
		
		SourceSDFWithSI(NumericType fluxDensityAJy, NumericType referenceFrequencyAHz, NumericType fluxDensityBJy, NumericType referenceFrequencyBHz) :
			/* Calculate the spectral index and flux density for 1 Hz frequency */
			_spectralIndex((log(fabs(fluxDensityAJy)) - log(fabs(fluxDensityBJy))) / (log(referenceFrequencyAHz) - log(referenceFrequencyBHz))),
		_refFreqA(referenceFrequencyAHz),
		_refFreqB(referenceFrequencyBHz)
		{
			_fluxDensityJy = fluxDensityAJy * std::pow((NumericType) 1.0 / referenceFrequencyAHz, _spectralIndex);
		}
		
		virtual NumericType FluxAtFrequency(NumericType frequencyHz) const
		{
			return _fluxDensityJy * std::pow(frequencyHz, _spectralIndex);
		}
		
		static NumericType FluxAtFrequency(NumericType fluxDensityAJy, NumericType referenceFrequencyAHz, NumericType fluxDensityBJy, NumericType referenceFrequencyBHz, NumericType requestedFrequency)
		{
			// if either fluxes are zero, or one of them is negative and the other not,
			// perform linear interpolation instead of power law interpolation
			bool signA = fluxDensityAJy < 0.0, signB = fluxDensityBJy < 0.0;
			if(fluxDensityAJy==0.0 || fluxDensityBJy==0.0 || (signA && !signB) || (signB && !signA))
			{
				NumericType slope =
					(fluxDensityBJy - fluxDensityAJy) /
					(referenceFrequencyBHz - referenceFrequencyAHz);
				return fluxDensityAJy + slope * (requestedFrequency - referenceFrequencyAHz);
			} else {
				NumericType si =
					log(fabs(fluxDensityBJy/fluxDensityAJy)) /
					log(referenceFrequencyBHz/referenceFrequencyAHz);
				return fluxDensityAJy * std::pow(requestedFrequency/referenceFrequencyAHz, si);
			}
		}
		
		virtual NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const
		{
			return _fluxDensityJy * (std::pow(endFrequency, _spectralIndex+1.0) - std::pow(startFrequency, _spectralIndex+1.0)) / ((_spectralIndex+1.0) * (endFrequency-startFrequency));
		}
		
		static NumericType IntegratedFlux(NumericType fluxDensityAJy, NumericType referenceFrequencyAHz, NumericType fluxDensityBJy, NumericType referenceFrequencyBHz, NumericType startFrequency, NumericType endFrequency)
		{
			// if either fluxes are zero, or one of them is negative and the other not,
			// perform linear interpolation instead of power law interpolation
			bool signA = fluxDensityAJy < 0.0, signB = fluxDensityBJy < 0.0;
			if(fluxDensityAJy==0.0 || fluxDensityBJy==0.0 || (signA && !signB) || (signB && !signA))
			{
				NumericType slope =
					(fluxDensityBJy - fluxDensityAJy) /
					(referenceFrequencyBHz - referenceFrequencyAHz);
				return slope * 0.5 * (
					endFrequency*endFrequency - startFrequency*startFrequency) /
					(endFrequency - startFrequency)
					+ (fluxDensityAJy - slope * referenceFrequencyAHz);
			} else {
				NumericType si =
					log(fluxDensityBJy/fluxDensityAJy) /
					log(referenceFrequencyBHz/referenceFrequencyAHz);
				return fluxDensityAJy * (
					std::pow(endFrequency/referenceFrequencyAHz, si) * endFrequency -
					std::pow(startFrequency/referenceFrequencyAHz, si) * startFrequency) /
					( (si+1.0) * (endFrequency - startFrequency) );
			}
		}
		
		NumericType SpectralIndex() const { return _spectralIndex; }
		
		virtual SourceSDF<NumericType> *Copy() const
		{
			SourceSDFWithSI<NumericType> *sdf = new SourceSDFWithSI<NumericType>();
			sdf->_fluxDensityJy = _fluxDensityJy;
			sdf->_spectralIndex = _spectralIndex;
			sdf->_refFreqA = _refFreqA;
			sdf->_refFreqB = _refFreqB;
			return sdf;
		}
		
		virtual std::string ToString() const
		{
			std::ostringstream s;
			s << "spectralindex "
				<< _refFreqA/1000000.0 << ' '
				<< FluxAtFrequency(_refFreqA) << ' ' << SpectralIndex();
			//  << FluxAtFrequency(_refFreqA) << ' ' << _refFreqA/1000000.0 << ' '
			//  << FluxAtFrequency(_refFreqB) << ' ' << _refFreqB/1000000.0;
			return s.str();
		}

	private:
		NumericType _fluxDensityJy, _spectralIndex;
		NumericType _refFreqA, _refFreqB;
};

#include "sourcesdfwithsamples.h"

template<typename NumericType>
inline SourceSDF<NumericType> *SourceSDF<NumericType>::ParseLine(std::vector<std::string>::const_iterator firstToken,
	std::vector<std::string>::const_iterator endToken)
{
	if(*firstToken == "spectralindex")
	{
		if(endToken - firstToken == 5)
		{
			NumericType
				fluxA = strtold((firstToken+1)->c_str(), 0),
				refFreqA = strtold((firstToken+2)->c_str(), 0)*1000000.0,
				fluxB = strtold((firstToken+3)->c_str(), 0),
				refFreqB = strtold((firstToken+4)->c_str(), 0)*1000000.0;
			return new SourceSDFWithSI<NumericType>(fluxA, refFreqA, fluxB, refFreqB);
		} else {
			NumericType
				refFreq = strtold((firstToken+1)->c_str(), 0)*1000000.0,
				flux = strtold((firstToken+2)->c_str(), 0),
				si = strtold((firstToken+3)->c_str(), 0);
			return new SourceSDFWithSI<NumericType>(flux, si, refFreq);
		}
	} else if(*firstToken == "sampled")
	{
		++firstToken;
		size_t fluxCount = atol(firstToken->c_str());
		SourceSDFWithSamples<NumericType> *sdf = new SourceSDFWithSamples<NumericType>();
		for(size_t i=0;i!=fluxCount;++i)
		{
			++firstToken;
			NumericType flux = strtold(firstToken->c_str(), 0);
			++firstToken;
			NumericType frequency = strtold(firstToken->c_str(), 0)*1000000.0;
			sdf->AddSample(flux, frequency);
		}
		return sdf;
	}
	else throw std::runtime_error("Could not parse sdf in line: invalid sdf type");
}

#endif
