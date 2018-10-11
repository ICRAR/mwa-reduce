#include "tilebeam2016.h"

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

#include <limits>

TileBeam2016::TileBeam2016(const double* delays, bool frequencyInterpolation, const std::string& searchPath) :
	Beam2016Implementation(delays, nullptr, searchPath),
	_frequencyInterpolation(frequencyInterpolation)
{
	for(size_t i=0;i!=16;++i)
	{
		_delaysInSteps[i] = delays[i];
	}
}

/**
	* Get the full Jones matrix response of the tile including the dipole
	* reponse and array factor incorporating any mutual coupling effects
	* from the impedance matrix. freq in Hz.
	*/
void TileBeam2016::getTabulatedResponse(double az, double za, double freq, std::complex< double >* result)
{
	// input are radians -> convert to degrees as implementation class expects :
	double az_deg = az*(180.00/M_PI);
	double za_deg = za*(180.00/M_PI);
	JonesMatrix jones = CalcJones( az_deg, za_deg, freq, 1 );
	result[0] = jones.j00;
	result[1] = jones.j01;
	result[2] = jones.j10;
	result[3] = jones.j11;
}

void TileBeam2016::getInterpolatedResponse(double az, double za, double freq, std::complex<double>* result)
{
	// input are radians -> convert to degrees as implementation class expects :
	double az_deg = az*(180.00/M_PI);
	double za_deg = za*(180.00/M_PI);

	JonesMatrix jones =  CalcJones( az_deg, za_deg, freq, 1 );
	result[0] = jones.j00;
	result[1] = jones.j01;
	result[2] = jones.j10;
	result[3] = jones.j11;
}
