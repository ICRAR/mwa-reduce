#include "tilebeam2016.h"

TileBeam2016::TileBeam2016(const double* delays, bool frequencyInterpolation, const std::string& searchPath) :
	Beam2016Implementation(delays, nullptr, searchPath),
	_frequencyInterpolation(frequencyInterpolation)
{ }

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
