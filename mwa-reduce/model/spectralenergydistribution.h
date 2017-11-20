#ifndef SPECTRAL_ENERGY_DISTRIBUTION_H
#define SPECTRAL_ENERGY_DISTRIBUTION_H

#include "measurement.h"

#include "../nlplfitter.h"
#include "../polarization.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <stdexcept>

//#define EXTRA_ASSERTIONS 1

class SpectralEnergyDistribution
{
public:
	virtual ~SpectralEnergyDistribution() { }
	virtual SpectralEnergyDistribution* Clone() const = 0;
	virtual std::string ToString() const = 0;
	virtual long double FluxAtFrequencyFromIndex(long double frequencyHz, size_t pIndex) const = 0;
	virtual long double IntegratedFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const = 0;
	virtual long double AverageFlux(long double startFrequency, long double endFrequency, PolarizationEnum polarization) const = 0;
	virtual bool operator<(const SpectralEnergyDistribution &other) const = 0;
	virtual void operator*=(double factor) = 0;
	virtual void operator+=(const SpectralEnergyDistribution &other) = 0;
	virtual long double ReferenceFrequencyHz() const = 0;
	
	long double FluxAtFrequency(long double frequencyHz, PolarizationEnum polarization) const
	{
#ifdef EXTRA_ASSERTIONS
		if(!Polarization::IsStokes(polarization))
			throw std::runtime_error("Cannot store specified polarization in model");
#endif
		return FluxAtFrequencyFromIndex(frequencyHz, Polarization::StokesToIndex(polarization));
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
			  log(std::abs(fluxDensityBJy/fluxDensityAJy)) /
				log(referenceFrequencyBHz/referenceFrequencyAHz);
			return fluxDensityAJy * std::pow(requestedFrequency/referenceFrequencyAHz, si);
		}
	}
	
	long double FluxAtChannel(size_t channelIndex, size_t channelCount, long double startFreq, long double endFreq, PolarizationEnum polarization) const
	{
		long double freq = startFreq + (long double) channelIndex * (endFreq - startFreq) / (long double) (channelCount-1);
		return FluxAtFrequency(freq, polarization);
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
		}
		else {
			long double si =
				log(fluxDensityBJy/fluxDensityAJy) /
				log(referenceFrequencyBHz/referenceFrequencyAHz);
			if(si == -1.0)
			{
				return (log(endFrequency) - log(startFrequency)) *
					(referenceFrequencyAHz * fluxDensityAJy) / (endFrequency - startFrequency);
			}
			else {
				return fluxDensityAJy * (
					std::pow(endFrequency/referenceFrequencyAHz, si) * endFrequency -
					std::pow(startFrequency/referenceFrequencyAHz, si) * startFrequency) /
					( (si+1.0) * (endFrequency - startFrequency) );
			}
		}
	}
	
};

#endif
