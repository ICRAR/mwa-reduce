#ifndef TILE_BEAM_H
#define TILE_BEAM_H

#include <complex>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>

namespace casa
{
	class MEpoch;
	class MPosition;
};

class TileBeam
{
public:
	TileBeam(const double *delays);
	
	void AnalyticGain(casa::MDirection &referenceDir, casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, double &x, double &y);
	
	void AnalyticGain(double raRad, double decRad, const casa::MDirection::Ref &ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double latitude, double frequencyHz, double &x, double &y);
	
	void AnalyticGain(double zenithAngle, double azimuth, double frequencyHz, double &x, double &y);
	
	void AnalyticJones(double raRad, double decRad, const casa::MDirection::Ref &ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double haZenith, double decZenith, double frequencyHz, std::complex<double>* gain);
	
	void AnalyticJones(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haZenith, double decZenith, std::complex<double> *gain);
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
