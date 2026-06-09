#ifndef TILE_BEAM_2016_H
#define TILE_BEAM_2016_H

#include <complex>
#include <map>
#include <set>

#include "beam2016implementation.h"

class TileBeam2016 : public Beam2016Implementation
{
public:
	TileBeam2016(const double *delays, bool frequencyInterpolation, const std::string& searchPath);
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, [[maybe_unused]] double ha, [[maybe_unused]] double dec, [[maybe_unused]] double haAntennaZenith, [[maybe_unused]] double decAntennaZenith, std::complex<double> *gain)
	{
		if(_frequencyInterpolation)
			getInterpolatedResponse(azimuth, zenithAngle, frequencyHz, gain);
		else
			getTabulatedResponse(azimuth, zenithAngle, frequencyHz, gain);
	}
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, std::complex<double> *gain)
	{
		if(_frequencyInterpolation)
			getInterpolatedResponse(azimuth, zenithAngle, frequencyHz, gain);
		else
			getTabulatedResponse(azimuth, zenithAngle, frequencyHz, gain);
	}
	#ifdef CUDA_SUPPORT
	void ArrayResponse(double *zenithAngle, double *azimuth, size_t npos, std::complex<double> *gain, size_t startChannel, size_t endChannel, size_t channelCount, size_t startFrequency, size_t endFrequency){
		getTabulatedResponse(azimuth, zenithAngle, npos, gain, startChannel, endChannel, channelCount, startFrequency, endFrequency);
	}
	#endif
private:
	bool _frequencyInterpolation;
	
	/**
	 * Get the full Jones matrix response of the tile including the dipole
	 * reponse and array factor incorporating any mutual coupling effects
	 * from the impedance matrix. freq in Hz.
	 */
	void getTabulatedResponse(double az, double za, double freq, std::complex<double>* result);
	
	#ifdef CUDA_SUPPORT
	void getTabulatedResponse(double *az, double *za, size_t npos, std::complex<double>* result, size_t startChannel, size_t endChannel, size_t channelCount, size_t startFrequency, size_t endFrequency);
	#endif
	/**
	 * Create a few tabulated responses and interpolated over these.
	 */
	void getInterpolatedResponse(double az, double za, double freq, std::complex<double>* result)
	{
		// Not implemented yet: just call normal function
		getTabulatedResponse(az, za, freq, result);
	}
};

#endif
