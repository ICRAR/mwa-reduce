#include "tilebeambase.h"
#include "tilebeam2013.h"
#include "tilebeam2014.h"

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/Measures/MeasConvert.h>

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_LATTITUDE = -26.703319; // Array latitude. degrees North

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_LONGITUDE = 116.67081;  // Array longitude. degrees East

template<typename Implementation>
const double TileBeamBase<Implementation>::MWA_HEIGHT = 377.0;         // Array altitude. meters above sea level

template<typename Implementation>
void TileBeamBase<Implementation>::ArrayResponse(casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, std::complex<double>* gain)
{
	casa::MeasFrame frame(arrayPos, time);
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MPosition wgs = casa::MPosition::Convert(arrayPos, casa::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casa::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
		
	casa::MDirection zenith(casa::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casa::MDirection zenithHaDec = casa::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getAngle().getValue()[0];
	double zenithDec = zenithHaDec.getAngle().getValue()[1];
	
	ArrayResponse(raRad, decRad, j2000Ref, j2000ToHaDec, j2000ToAzelGeo, arrLatitude, zenithHa, zenithDec, frequencyHz, gain);
}

template<typename Implementation>
void TileBeamBase<Implementation>::ArrayResponse(double raRad, double decRad, const casa::MDirection::Ref &j2000Ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double haAntennaZenith, double decAntennaZenith, double frequencyHz, std::complex<double>* gain)
{
	static const casa::Unit radUnit("rad");
	casa::MDirection imageDir(casa::MVDirection(
		casa::Quantity(raRad, radUnit),     // RA
		casa::Quantity(decRad,radUnit)),  // DEC
		j2000Ref);
	
	// convert ra, dec to ha
	casa::MDirection hadec = j2000ToHaDec(imageDir);
	double ha = hadec.getValue().get()[0];
	double sinLat, cosLat;
	sincos(arrLatitude, &sinLat, &cosLat);
	double sinDec, cosDec;
	sincos(decRad, &sinDec, &cosDec);
	double cosHA = cos(ha);
	double zenithDistance = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
	casa::MDirection azel = j2000ToAzelGeo(imageDir);
	double azimuth = azel.getValue().get()[0];
	//std::cout << "ha=" << ha*180.0/M_PI << '-' << haAntennaZenith*180.0/M_PI << ", za=" << zenithDistance*180.0/M_PI << ", az=" << azimuth*180.0/M_PI << '\n';
	
	ArrayResponse(zenithDistance, azimuth, frequencyHz, ha, decRad, haAntennaZenith, decAntennaZenith, gain);
}

template<typename Implementation>
void TileBeamBase<Implementation>::PrecalculatePositionInfo(TileBeamBase::PrecalcPosInfo& posInfo, casa::MEpoch& time, casa::MPosition& arrayPos, double raRad, double decRad)
{
	casa::MeasFrame frame(arrayPos, time);
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MPosition wgs = casa::MPosition::Convert(arrayPos, casa::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casa::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
		
	casa::MDirection zenith(casa::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casa::MDirection zenithHaDec = casa::MDirection::Convert(zenith, hadecRef)();
	posInfo.haAntennaZenith = zenithHaDec.getAngle().getValue()[0];
	posInfo.decAntennaZenith = zenithHaDec.getAngle().getValue()[1];
	
	casa::MDirection imageDir(casa::MVDirection(raRad, decRad), j2000Ref);
	
	// convert ra, dec to ha
	casa::MDirection hadec = j2000ToHaDec(imageDir);
	posInfo.ha = hadec.getValue().get()[0];
	posInfo.dec = decRad;
	double sinLat, cosLat;
	sincos(arrLatitude, &sinLat, &cosLat);
	double sinDec, cosDec;
	sincos(decRad, &sinDec, &cosDec);
	double cosHA = cos(posInfo.ha);
	posInfo.zenithAngle = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
	casa::MDirection azel = j2000ToAzelGeo(imageDir);
	posInfo.azimuth = azel.getValue().get()[0];
}

template class TileBeamBase<TileBeam2013>;

template class TileBeamBase<TileBeam2014>;

