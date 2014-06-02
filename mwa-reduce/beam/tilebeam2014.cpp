#include "tilebeam2014.h"
#include "joneslookupdipole.h"

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

TileBeam2014::TileBeam2014(const double* delays) :
	_dipoleHeight(0.28), /* Seems to be 0.3 in the RTS, 0.278 in beam script */
	_dipoleSeparations(1.100),
	_delayStep(435.0e-12)
{
	const double dipoleEast[16] = {-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5};
	const double dipoleNorth[16] = {1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5};
	for(size_t i=0;i!=16;++i)
	{
		_dipoleNorth[i] = dipoleNorth[i] * _dipoleSeparations;
		_dipoleEast[i] = dipoleEast[i] * _dipoleSeparations;
		_delays[i] = delays[i]*SPEED_OF_LIGHT*_delayStep;
	}
	JonesLookupDipole::Initialize();
	TileImpedance::Initialize();
}

void TileBeam2014::invert32x32(const std::complex<double>* input, std::complex<double>* output)
{
	gsl_matrix_complex
		*gslInput = gsl_matrix_complex_alloc(32, 32),
		*gslOutput = gsl_matrix_complex_alloc(32, 32);
	gsl_permutation *perm = gsl_permutation_alloc(32);
	
	memcpy(gsl_matrix_complex_ptr(gslInput, 0, 0), input, sizeof(std::complex<double>)*32*32);
	
		// Make LU decomposition of matrix m
	int s;
	gsl_linalg_complex_LU_decomp(gslInput, perm, &s);

	// Invert the matrix m
	gsl_linalg_complex_LU_invert(gslInput, perm, gslOutput);
	
	memcpy(output, gsl_matrix_complex_ptr(gslOutput, 0, 0), sizeof(std::complex<double>)*32*32);
	
	gsl_matrix_complex_free(gslInput);
	gsl_matrix_complex_free(gslOutput);
	gsl_permutation_free(perm);
}

/**
	* Get the full Jones matrix response of the tile including the dipole
	* reponse and array factor incorporating any mutual coupling effects
	* from the impedance matrix. freq in Hz.
	*/
void TileBeam2014::getResponse(double az, double za, double freq, std::complex< double >* result)
{
	const FrequencyCacheInfo *cacheInfo;
	if(_frequencyCache.count(freq) == 0)
	{
		// Frequency is not in cache: fill cache
		FrequencyCacheInfo newCacheInfo;
		getPortCurrents(freq, newCacheInfo.current, _delays);
		newCacheInfo.lambda = SPEED_OF_LIGHT / freq;
		double zenithDelays[16];
		for(size_t i=0; i!=16; ++i) zenithDelays[i] = 0.0;
		getArrayFactor(0.0, 0.0, freq, newCacheInfo.zax, newCacheInfo.zay, zenithDelays);
		
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

	const bool useLookup = true;
	if(useLookup)
	{
		std::complex<double> dipoleJones[4];
		JonesLookupDipole::Interpolate(dipoleJones, az, za, freq);
		result[0] = ax * dipoleJones[0];
		result[1] = ax * dipoleJones[1];
		result[2] = ay * dipoleJones[2];
		result[3] = ay * dipoleJones[3];
	}
	else {
		double dipoleJones[4];
		getJonesShortDipole(az, za, freq, dipoleJones);
		result[0] = ax * dipoleJones[0];
		result[1] = ax * dipoleJones[1];
		result[2] = ay * dipoleJones[2];
		result[3] = ay * dipoleJones[3];
	}
}
