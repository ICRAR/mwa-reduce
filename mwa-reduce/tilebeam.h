#ifndef TILE_BEAM_H
#define TILE_BEAM_H

#include <complex>

class TileBeam
{
public:
	TileBeam(const double *delays);
		
	void AnalyticGain(double zenithAngle, double azimuth, double frequencyHz, double &x, double &y);
private:
	double _dipoleSize; // height of dipole
	double _dipoleSeparations;
	double _delayStep;
	bool _zenithNorm;
	double _zenithNormFactor;
	
	double _dipoleNorth[16];
	double _dipoleEast[16];
	double _dipoleHeight[16];
	double _delays[16];
};

#endif
