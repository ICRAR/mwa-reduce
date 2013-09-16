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
	struct PrecalcPosInfo
	{
		double zenithAngle;
		double azimuth;
		double ha;
		double dec;
		double haAntennaZenith;
		double decAntennaZenith;
	};
	
	TileBeam(const double *delays);
	
	void AnalyticGain(casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, double &x, double &y);
	
	void AnalyticGain(double raRad, double decRad, const casa::MDirection::Ref &j2000Ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double latitude, double frequencyHz, double &x, double &y);
	
	void AnalyticGain(double zenithAngle, double azimuth, double frequencyHz, double &x, double &y);
	
	void AnalyticJones(casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, std::complex<double>* gain);
	
	void AnalyticJones(double raRad, double decRad, const casa::MDirection::Ref &j2000Ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double haZenith, double decZenith, double frequencyHz, std::complex<double>* gain);
	
	void AnalyticJones(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain);
	
	void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad);
	
	void AnalyticJones(const PrecalcPosInfo& posInfo, double frequencyHz, std::complex<double> *gain)
	{
		AnalyticJones(
			posInfo.zenithAngle,
			posInfo.azimuth,
			frequencyHz,
			posInfo.ha,
			posInfo.dec,
			posInfo.haAntennaZenith,
			posInfo.decAntennaZenith,
			gain
		);
	}
	
	const static double MWA_LATTITUDE; // Array latitude. degrees North
	const static double MWA_LONGITUDE; // Array longitude. degrees East
	const static double MWA_HEIGHT;    // Array altitude. meters above sea level
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
