#ifndef TILE_BEAM_2014_H
#define TILE_BEAM_2014_H

#include <complex>
#include <map>
#include <set>

#include "tileimpedance.h"
#include "lnaimpedance.h"

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s
#endif

class TileBeam2014
{
public:
	TileBeam2014(const double *delays, bool frequencyInterpolation, const std::string& searchPath);
	
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
	
private:
	double _dipoleEast[16];
	double _dipoleNorth[16];
	double _delays[16];
	const double _dipoleHeight, _dipoleSeparations, _delayStep;
	std::set<double> _tabulationFrequencies;
	bool _frequencyInterpolation;
	static bool _firstInit;
	
	struct FrequencyCacheInfo
	{
		std::complex<double> current[32];
		std::complex<double> zax, zay;
		double lambda;
	};
	
	/**
		* Calculate the Jones matrix for a short dipole.
		* This is defined by purely geometric projection of unit vectors
		* on the sky onto the unit vector defined by the dipole's direction.
	*/
	void getJonesShortDipole(double az, double za, double freq, double* result)
	{
		// apply the groundscreen factor, which is independent of az
		double lambda = SPEED_OF_LIGHT / freq;
		double cosZa = cos(za), sinAz, cosAz;
		sincos(az, &sinAz, &cosAz);
		double gs = groundScreen(cosZa, lambda);
		gs /= groundScreenZenith(lambda);
		result[0] = cosZa*sinAz*gs;
		result[1] = cosAz*gs;
		result[2] = cosZa*cosAz*gs;
		result[3] = -sinAz*gs;
	}
	
	/**
	 * Calculate the groundscreen effect for an ideal infinite groundscreen
	 * given the dipole's height above the screen and the frequency (Hz)
	 */
	double groundScreen(double cosZA, double lambda) const
	{
		return sin(M_PI * (2.0*_dipoleHeight/lambda) * cosZA)*2.0;
	}
	
	/**
	 * Calculate the groundscreen effect for an ideal infinite groundscreen
	 * given the dipole's height above the screen and the frequency (Hz)
	 */
	double groundScreenZenith(double lambda) const
	{
		return sin(M_PI * (2.0*_dipoleHeight/lambda))*2.0;
	}
	
	/**
	 * Return the port currents on a tile given the freq (Hz) and delays
	 */
 	static void getPortCurrents(double freq, std::complex<double>* current, const double* delays);
	
	/**
	 * Get the scalar array factor response of the array for a given
	 * freq (Hz) and delay settings.
	 * az and za (radian) are numpy arrays of equal length defining a set
	 * of points to calculate the response for.
	 * delays is a 2D array of integer delay steps for the Y and X pol
	 * respectively.
	 * 
	 * Result are in same coords as the az/za input arrays
	 */
  void getArrayFactor(double az, double za, double freq, std::complex<double>& ax, std::complex<double>& ay, const double* delays)
	{
		double lambda = SPEED_OF_LIGHT / freq;
		std::complex<double> portCurrent[32];
		getPortCurrents(freq, portCurrent, delays);
		
		// now calculate the array factor using these port currents
		double sz = sin(za);
		double kx = (2.0*M_PI/lambda)*sin(az)*sz;
		double ky = (2.0*M_PI/lambda)*cos(az)*sz;
		
		ax = 0.0; ay = 0.0;
		
		// Only calculate if above horizon
		if(za < M_PI/2.0)
		{
			for(size_t i=0; i!=16; ++i)
			{
				double ph = kx*_dipoleEast[i] + ky*_dipoleNorth[i];
				double s, c;
				sincos(ph, &s, &c);
				ax += portCurrent[i + 16] * std::complex<double>(c, s); // X dipoles
				ay += portCurrent[     i] * std::complex<double>(c, s); // Y dipoles
			}
		}
	}

	/**
	 * Get the scalar array factor response of the array for a given
	 * freq (Hz) and delay settings.
	 */
  void getArrayFactor(double az, double za, const FrequencyCacheInfo& cacheInfo, std::complex<double>& ax, std::complex<double>& ay)
	{
		// now calculate the array factor using these port currents
		double sz = sin(za);
		double sinAz, cosAz;
		sincos(az, &sinAz, &cosAz);
		double kx = (2.0*M_PI/cacheInfo.lambda) * sinAz * sz;
		double ky = (2.0*M_PI/cacheInfo.lambda) * cosAz * sz;
		
		ax = 0.0; ay = 0.0;
		
		// Only calculate if above horizon
		if(za < M_PI/2.0)
		{
			for(size_t i=0; i!=16; ++i)
			{
				double ph = kx*_dipoleEast[i] + ky*_dipoleNorth[i];
				double s, c;
				sincos(ph, &s, &c);
				ax += cacheInfo.current[i + 16] * std::complex<double>(c, s); // X dipoles
				ay += cacheInfo.current[     i] * std::complex<double>(c, s); // Y dipoles
			}
		}
	}

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
        
	static void invert32x32(const std::complex<double>* input, std::complex<double>* output);
	
	std::map<double, FrequencyCacheInfo> _frequencyCache;
};

#undef SPEED_OF_LIGHT

#endif
