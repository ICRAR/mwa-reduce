#include "tilebeam.h"

#include <cmath>
#include <cstring>

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

// Based on code from Daniel Mitchel
// 2012-02-13
// taken from the RTS codebase
// Optimized 2012-11-17 by Offringa.

TileBeam::TileBeam(const double *delays) :
	_dipoleSize(0.278),
	_dipoleSeparations(1.1),
	_delayStep(435.0e-12),
	_zenithNorm(true)
{
	const double dipoleNorth[16] = {1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5};
	const double dipoleEast[16] = {-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5};
	const double dipoleHeight[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	for(size_t i=0;i!=16;++i)
	{
		_dipoleNorth[i] = dipoleNorth[i] * _dipoleSeparations;
		_dipoleEast[i] = dipoleEast[i] * _dipoleSeparations;
		_dipoleHeight[i] = dipoleHeight[i] * _dipoleSeparations;
		_delays[i] = delays[i]*SPEED_OF_LIGHT*_delayStep;
	}
}

void TileBeam::AnalyticGain(double zenithAngle, double azimuth, double frequencyHz, std::complex<double> &x, std::complex<double> &y)
{
	// direction cosines (relative to zenith) for direction az,za
	const double sinZenith = sin(zenithAngle), cosZenith = cos(zenithAngle);
	const double projectionEast = sinZenith * sin(azimuth);
	const double projectionNorth = sinZenith * cos(azimuth);
	const double projectionHeight = cosZenith;
	
	const double lambda = SPEED_OF_LIGHT / frequencyHz;
	const double twoPiOverLambda = 2.0 * M_PI / lambda;

	// dipole position within the tile
	std::complex<double> arrayFactor = 0.0;
	for(size_t i=0;i!=16;++i)
	{
		// relative dipole phase for a source at (theta,phi)
		double rotation = twoPiOverLambda*(_dipoleEast[i]*projectionEast + _dipoleNorth[i]*projectionNorth +
			_dipoleHeight[i]*projectionHeight - _delays[i]);
		double rotSin, rotCos;
		sincos(rotation, &rotSin, &rotCos);
    arrayFactor += std::complex<double>(rotCos, rotSin);
	}
	arrayFactor /= 16.0;

	double groundPlane;
	
  // make sure we filter out the bottom hemisphere
	if(zenithAngle > M_PI)
		groundPlane = 0.0;
	else
		groundPlane = 2.0 * sin(twoPiOverLambda * _dipoleSize * cosZenith);
	
	// normalize to zenith
	if(_zenithNorm)
		groundPlane /= 2.0 * sin(twoPiOverLambda * _dipoleSize);

	// response of the 2 tile polarizations
	// gains due to forshortening
	double dipole_ns = sqrt(1.0 - projectionNorth*projectionNorth);
	double dipole_ew = sqrt(1.0 - projectionEast*projectionEast);

	// voltage responses of the polarizations from an unpolarized source
	// this is effectively the YY voltage gain
	y = dipole_ns * groundPlane * arrayFactor; // gain_ns
	// this is effectively the XX voltage gain
	x = dipole_ew * groundPlane * arrayFactor; // gain_ew
	
}
