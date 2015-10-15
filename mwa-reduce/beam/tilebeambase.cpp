#include "tilebeambase.h"
#include "tilebeam2013.h"
#include "tilebeam2014.h"

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MeasConvert.h>

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_LATTITUDE = -26.703319; // Array latitude. degrees North

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_LONGITUDE = 116.67081;  // Array longitude. degrees East

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_HEIGHT = 377.0;         // Array altitude. meters above sea level

template<typename Implementation>
void TileBeamBase<Implementation>::ArrayResponse(casacore::MEpoch &time, casacore::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, std::complex<double>* gain)
{
	casacore::MeasFrame frame(arrayPos, time);
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, frame);
	const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MPosition wgs = casacore::MPosition::Convert(arrayPos, casacore::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casacore::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
		
	casacore::MDirection zenith(casacore::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casacore::MDirection zenithHaDec = casacore::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getAngle().getValue()[0];
	double zenithDec = zenithHaDec.getAngle().getValue()[1];
	
	ArrayResponse(raRad, decRad, j2000Ref, j2000ToHaDec, j2000ToAzelGeo, arrLatitude, zenithHa, zenithDec, frequencyHz, gain);
}

template<typename Implementation>
void TileBeamBase<Implementation>::ArrayResponse(double raRad, double decRad, const casacore::MDirection::Ref &j2000Ref, casacore::MDirection::Convert &j2000ToHaDec, casacore::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double haAntennaZenith, double decAntennaZenith, double frequencyHz, std::complex<double>* gain)
{
	static const casacore::Unit radUnit("rad");
	casacore::MDirection imageDir(casacore::MVDirection(
		casacore::Quantity(raRad, radUnit),     // RA
		casacore::Quantity(decRad,radUnit)),  // DEC
		j2000Ref);
	
	// convert ra, dec to ha
	casacore::MDirection hadec = j2000ToHaDec(imageDir);
	double ha = hadec.getValue().get()[0];
	double sinLat, cosLat;
	sincos(arrLatitude, &sinLat, &cosLat);
	double sinDec, cosDec;
	sincos(decRad, &sinDec, &cosDec);
	double cosHA = cos(ha);
	double zenithDistance = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
	casacore::MDirection azel = j2000ToAzelGeo(imageDir);
	double azimuth = azel.getValue().get()[0];
	//std::cout << "ha=" << ha*180.0/M_PI << '-' << haAntennaZenith*180.0/M_PI << ", za=" << zenithDistance*180.0/M_PI << ", az=" << azimuth*180.0/M_PI << '\n';
	
	ArrayResponse(zenithDistance, azimuth, frequencyHz, ha, decRad, haAntennaZenith, decAntennaZenith, gain);
}

template<typename Implementation>
void TileBeamBase<Implementation>::PrecalculatePositionInfo(TileBeamBase::PrecalcPosInfo& posInfo, casacore::MEpoch& time, casacore::MPosition& arrayPos, double raRad, double decRad)
{
	casacore::MeasFrame frame(arrayPos, time);
	const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
	const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO, frame);
	const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MPosition wgs = casacore::MPosition::Convert(arrayPos, casacore::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casacore::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
		
	casacore::MDirection zenith(casacore::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casacore::MDirection zenithHaDec = casacore::MDirection::Convert(zenith, hadecRef)();
	posInfo.haAntennaZenith = zenithHaDec.getAngle().getValue()[0];
	posInfo.decAntennaZenith = zenithHaDec.getAngle().getValue()[1];
	
	casacore::MDirection imageDir(casacore::MVDirection(raRad, decRad), j2000Ref);
	
	// convert ra, dec to ha
	casacore::MDirection hadec = j2000ToHaDec(imageDir);
	posInfo.ha = hadec.getValue().get()[0];
	posInfo.dec = decRad;
	double sinLat, cosLat;
	sincos(arrLatitude, &sinLat, &cosLat);
	double sinDec, cosDec;
	sincos(decRad, &sinDec, &cosDec);
	double cosHA = cos(posInfo.ha);
	posInfo.zenithAngle = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
	casacore::MDirection azel = j2000ToAzelGeo(imageDir);
	posInfo.azimuth = azel.getValue().get()[0];
}

template class TileBeamBase<TileBeam2013>;

template class TileBeamBase<TileBeam2014>;

