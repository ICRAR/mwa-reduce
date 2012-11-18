#ifndef SOURCE_STRENGTH_H
#define SOURCE_STRENGTH_H

#include <cmath>

template<typename NumericType=long double>
class SourceStrength
{
	public:
		SourceStrength() :
		_fluxDensityJy(0.0),
		_spectralIndex(0.0)
		{
		}
		
		SourceStrength(NumericType fluxDensityJy, NumericType spectralIndex, NumericType siReferenceFrequencyHz) :
		/* Calculate the flux density for 1 Hz frequency */
			_fluxDensityJy( fluxDensityJy *
				std::pow((NumericType) 1.0 / siReferenceFrequencyHz, spectralIndex) ),
			_spectralIndex( spectralIndex )
		{
		}
		
		SourceStrength(NumericType fluxDensityAJy, NumericType referenceFrequencyAHz, NumericType fluxDensityBJy, NumericType referenceFrequencyBHz) :
			/* Calculate the spectral index and flux density for 1 Hz frequency */
		_spectralIndex( (log(fabs(fluxDensityAJy)) - log(fabs(fluxDensityBJy))) / (log(referenceFrequencyAHz) - log(referenceFrequencyBHz)) )
		{
			_fluxDensityJy = fluxDensityAJy * std::pow((NumericType) 1.0 / referenceFrequencyAHz, _spectralIndex);
		}
		
		NumericType FluxAtFrequency(NumericType frequencyHz) const
		{
			return _fluxDensityJy * std::pow(frequencyHz, _spectralIndex);
		}
		
		NumericType FluxAtFrequency(size_t channelIndex, size_t channelCount, NumericType startFreq, NumericType endFreq) const
		{
			NumericType freq = startFreq + NumericType(channelIndex) * (endFreq - startFreq) / NumericType(channelCount-1);
			return _fluxDensityJy * std::pow(freq, _spectralIndex);
		}
		
		NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const
		{
			return _fluxDensityJy * (std::pow(endFrequency, _spectralIndex+1.0) - std::pow(startFrequency, _spectralIndex+1.0)) / ((_spectralIndex+1.0) * (endFrequency-startFrequency));
		}
		NumericType SpectralIndex() const { return _spectralIndex; }
	private:
		NumericType _fluxDensityJy, _spectralIndex;
};

#endif
