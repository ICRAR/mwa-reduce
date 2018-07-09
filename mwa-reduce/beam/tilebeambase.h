#ifndef TILE_BEAM_BASE_H
#define TILE_BEAM_BASE_H

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>

namespace casacore
{
	class MEpoch;
	class MPosition;
};

template<typename Implementation>
class TileBeamBase : private Implementation
{
public:
	TileBeamBase(const double *delays, bool frequencyInterpolation = true) : Implementation(delays, frequencyInterpolation) { }
	
	struct PrecalcPosInfo
	{
		double zenithAngle;
		double azimuth;
		double ha;
		double dec;
		double haAntennaZenith;
		double decAntennaZenith;
	};
	
	void ArrayResponse(casacore::MEpoch &time, casacore::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, std::complex<double>* gain);
	
	void ArrayResponse(double raRad, double decRad, const casacore::MDirection::Ref &j2000Ref, casacore::MDirection::Convert &j2000ToHaDec, casacore::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double haZenith, double decZenith, double frequencyHz, std::complex<double>* gain);
	
	void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, casacore::MEpoch &time, casacore::MPosition &arrayPos, double raRad, double decRad);
	
	void ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain)
	{
		Implementation::ArrayResponse(zenithAngle, azimuth, frequencyHz, ha, dec, haAntennaZenith, decAntennaZenith, gain);
	}
	
	void ArrayResponse(const PrecalcPosInfo& posInfo, double frequencyHz, std::complex<double> *gain)
	{
		ArrayResponse(
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
	
private:
	const static double MWA_LATTITUDE; // Array latitude. degrees North
	const static double MWA_LONGITUDE; // Array longitude. degrees East
	const static double MWA_HEIGHT;    // Array altitude. meters above sea level
	
};

#endif
