#ifndef TILE_BEAM_2014_H
#define TILE_BEAM_2014_H

#include <complex>
#include <map>

#include "tileimpedance.h"
#include "lnaimpedance.h"

#ifndef SPEED_OF_LIGHT
#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s
#endif

class TileBeam2014
{
public:
	TileBeam2014(const double *delays);
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain)
	{
		getResponse(azimuth, zenithAngle, frequencyHz, gain);
	}
	
private:
	double _dipoleEast[16];
	double _dipoleNorth[16];
	double _delays[16];
	const double _dipoleHeight, _dipoleSeparations, _delayStep;
	
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
		double gs = groundScreen(za, freq);
		gs /= groundScreen(0.0, freq);
		result[0] = cos(za)*sin(az)*gs;
		result[1] = -cos(az)*gs;
		result[2] = cos(za)*cos(az)*gs;
		result[3] = -sin(az)*gs;
	}
	
	/**
	 * Calculate the groundscreen effect for an ideal infinite groundscreen
	 * given the dipole's height above the screen and the frequency (Hz)
	 */
	double groundScreen(double za, double lambda) const
	{
		return sin(M_PI * (2.0*_dipoleHeight/lambda) * cos(za))*2.0;
	}
	
	/**
	 * Return the port currents on a tile given the freq (Hz) and delays (integer)
	 */
 	void getPortCurrents(double freq, std::complex<double>* current)
	{
		double lambda = SPEED_OF_LIGHT / freq;
		std::complex<double> ph_rot[32];
		for(size_t i=0; i!=16; ++i) {
			double phase = -2.0 * M_PI * _delays[i] / lambda;
			double s, c;
			sincos(phase, &s, &c);
			ph_rot[i] = std::complex<double>(c, s);
			ph_rot[i+16] = ph_rot[i];
		}
		
		std::complex<double> lnaImp = LNAImpedance::Get(freq);
		std::complex<double> zTotal[32*32];
		TileImpedance::Get(freq, zTotal);
		// Add lna impedance to diagonal values
		for(size_t diag=0; diag!=32; ++diag)
			zTotal[diag*33] += lnaImp;
		std::complex<double> invertedZ[32*32], *invertedZPtr = invertedZ;
		invert32x32(zTotal, invertedZ);
		// Compute inner product invertedZ . ph_rot
		for(size_t j=0; j!=32; ++j)
		{
			for(size_t i=0; i!=32; ++i) {
				current[j] += *invertedZPtr * ph_rot[i];
				++invertedZPtr;
			}
		}
	}
	
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
  void getArrayFactor(double az, double za, double freq, std::complex<double>& ax, std::complex<double>& ay)
	{
		double lambda = SPEED_OF_LIGHT / freq;
		std::complex<double> portCurrent[32];
		getPortCurrents(freq, portCurrent);
		
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
	 * az and za (radian) are numpy arrays of equal length defining a set
	 * of points to calculate the response for.
	 * delays is a 2D array of integer delay steps for the Y and X pol
	 * respectively.
	 * 
	 * Result are in same coords as the az/za input arrays
	 */
  void getArrayFactor(double az, double za, const FrequencyCacheInfo& cacheInfo, std::complex<double>& ax, std::complex<double>& ay)
	{
		// now calculate the array factor using these port currents
		double sz = sin(za);
		double kx = (2.0*M_PI/cacheInfo.lambda)*sin(az)*sz;
		double ky = (2.0*M_PI/cacheInfo.lambda)*cos(az)*sz;
		
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
	 * from the impedance matrix. freq in Hz, delays in unit steps of
	 * beamformer delays.
	 */
	void getResponse(double az, double za, double freq, std::complex<double>* result)
	{
		const FrequencyCacheInfo *cacheInfo;
		if(_frequencyCache.count(freq) == 0)
		{
			// Frequency is not in cache: fill cache
			FrequencyCacheInfo newCacheInfo;
			getArrayFactor(0.0, 0.0, freq, newCacheInfo.zax, newCacheInfo.zay);			
			getPortCurrents(freq, newCacheInfo.current);
			newCacheInfo.lambda = SPEED_OF_LIGHT / freq;
			
			_frequencyCache.insert(std::make_pair(freq, newCacheInfo));
			cacheInfo = &newCacheInfo;
		}
		else {
			cacheInfo = &_frequencyCache.find(freq)->second;
		}
			
		std::complex<double> ax, ay;
		getArrayFactor(az, za, *cacheInfo, ax, ay);
		
		// zenith response to normalise to
		ax /= std::abs(cacheInfo->zax);
		ay /= std::abs(cacheInfo->zay);

		double dipoleJones[4];
		getJonesShortDipole(az, za, freq, dipoleJones);
		result[0] = ax * dipoleJones[0];
		result[1] = ax * dipoleJones[1];
		result[2] = ay * dipoleJones[2];
		result[3] = ay * dipoleJones[3];
	}
        
	void invert32x32(const std::complex<double>* input, std::complex<double>* output);
	
	std::map<double, FrequencyCacheInfo> _frequencyCache;
};

#undef SPEED_OF_LIGHT

#endif
