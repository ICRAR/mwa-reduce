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
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain)
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
	
private:
	double _delaysInSteps[16];
	bool _frequencyInterpolation;
	
	/**
	 * Get the full Jones matrix response of the tile including the dipole
	 * reponse and array factor incorporating any mutual coupling effects
	 * from the impedance matrix. freq in Hz.
	 */
	void getTabulatedResponse(double az, double za, double freq, std::complex<double>* result);
	
	/**
	 * Create a few tabulated responses and interpolated over these.
	 */
	void getInterpolatedResponse(double az, double za, double freq, std::complex<double>* result);
};

#endif
