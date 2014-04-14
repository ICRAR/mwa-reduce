#ifndef TILE_BEAM_2013_H
#define TILE_BEAM_2013_H

#include <complex>

class TileBeam2013
{
public:
	TileBeam2013(const double *delays);
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain);
	
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
