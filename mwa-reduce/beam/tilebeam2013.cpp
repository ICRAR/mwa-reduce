#include "tilebeam.h"

#include <cmath>
#include <cstring>

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

// Based on code from Daniel Mitchel
// 2012-02-13
// taken from the RTS codebase
// Optimized 2012-11-17 by Offringa.

TileBeam2013::TileBeam2013(const double *delays) :
	_dipoleSize(0.278), /* Seems to be 0.3 in the RTS, 0.278 in beam script */
	_dipoleSeparations(1.100),
	_delayStep(435.0e-12),
	_zenithNorm(true)
{
	const double dipoleNorth[16] = {1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5};
	const double dipoleEast[16] = {-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5};
	//const double dipoleEast[16] = {-1.5,-1.5,-1.5,-1.5,-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,1.5,1.5,1.5,1.5};
	//const double dipoleNorth[16] = {-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5};
	const double dipoleHeight[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	for(size_t i=0;i!=16;++i)
	{
		_dipoleNorth[i] = dipoleNorth[i] * _dipoleSeparations;
		_dipoleEast[i] = dipoleEast[i] * _dipoleSeparations;
		_dipoleHeight[i] = dipoleHeight[i] * _dipoleSeparations;
		_delays[i] = delays[i]*SPEED_OF_LIGHT*_delayStep;
	}
}

/*
void TileBeam::AnalyticGain(casa::MEpoch &time, casa::MPosition &arrayPos, double raRad, double decRad, double frequencyHz, double &x, double &y)
{
	casa::MeasFrame frame(arrayPos, time);
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MPosition wgs = casa::MPosition::Convert(arrayPos, casa::MPosition::WGS84)();
	double latitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casa::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
		
	AnalyticGain(raRad, decRad, j2000Ref, j2000ToHaDec, j2000ToAzelGeo, latitude, frequencyHz, x, y);
}

void TileBeam::AnalyticGain(double raRad, double decRad, const casa::MDirection::Ref &ref, casa::MDirection::Convert &j2000ToHaDec, casa::MDirection::Convert &j2000ToAzelGeo, double arrLatitude, double frequencyHz, double &x, double &y)
{
	static const casa::Unit radUnit("rad");
	casa::MDirection imageDir(casa::MVDirection(
		casa::Quantity(raRad, radUnit),     // RA
		casa::Quantity(decRad,radUnit)),  // DEC
		ref);
	
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
	
	AnalyticGain(zenithDistance, azimuth, frequencyHz, x, y);
}

void TileBeam::AnalyticGain(double zenithAngle, double azimuth, double frequencyHz, double &x, double &y)
{
	// direction cosines (relative to zenith) for direction az,za
	double sinZenith, cosZenith, sinAzimuth, cosAzimuth;
	sincos(zenithAngle, &sinZenith, &cosZenith);
	sincos(azimuth, &sinAzimuth, &cosAzimuth);

	const double projectionEast = sinZenith * sinAzimuth;
	const double projectionNorth = sinZenith * cosAzimuth;
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
	
	double arrPower = (arrayFactor.real()*arrayFactor.real() + arrayFactor.imag()*arrayFactor.imag()) * groundPlane;
	y = dipole_ns * arrPower; // gain_ns
	// this is effectively the XX voltage gain
	x = dipole_ew * arrPower; // gain_ew
	
}*/

void TileBeam2013::ArrayResponse(double zenithAngle, double azimuth, double frequencyHz, double ha, double dec, double haAntennaZenith, double decAntennaZenith, std::complex<double> *gain)
{
	// direction cosines (relative to zenith) for direction az,za
	double sinZenith, cosZenith, sinAzimuth, cosAzimuth;
	sincos(zenithAngle, &sinZenith, &cosZenith);
	sincos(azimuth, &sinAzimuth, &cosAzimuth);

	const double projectionEast = sinZenith * sinAzimuth;
	const double projectionNorth = sinZenith * cosAzimuth;
	const double projectionHeight = cosZenith;
	
	const double lambda = SPEED_OF_LIGHT / frequencyHz;
	const double twoPiOverLambda = 2.0 * M_PI / lambda;
	
	/*const double cPhi = cosAzimuth, sPhi = sinAzimuth;
	const double cTheta = cosZenith, sTheta = -sinZenith;
	const double cPsi = cosAzimuth, sPsi = -sinAzimuth;
	
	const double factEE = cPsi*cPhi - cTheta*sPhi*sPsi;
	const double factEN = cPsi*sPhi + cTheta*cPhi*sPsi;
	const double factEH = sPsi*sTheta;
	
	const double factNE = -sPsi*cPhi - cTheta*sPhi*cPsi;
	const double factNN = -sPsi*sPhi + cTheta*cPhi*cPsi;
	const double factNH = cPsi*sTheta;
	
	const double factHE = sTheta*sPhi;
	const double factHN = -sTheta*cPhi;
	const double factHH = cTheta;*/

	// dipole position within the tile
	std::complex<double> arrayFactor = 0.0;
	for(size_t i=0;i!=16;++i)
	{
		// relative dipole phase for a source at (theta,phi)
		double rotation = twoPiOverLambda*(_dipoleEast[i]*projectionEast + _dipoleNorth[i]*projectionNorth +
			_dipoleHeight[i]*projectionHeight - _delays[i]);
		double rotSin, rotCos;
		sincos(rotation, &rotSin, &rotCos);
		std::complex<double> phaseShift = std::complex<double>(rotCos, rotSin);
    arrayFactor += phaseShift;
	}
	arrayFactor /= 16.0;

	double groundPlane;
	
  // make sure we filter out the bottom hemisphere
	if(zenithAngle > M_PI*0.5 || zenithAngle < -M_PI*0.5)
		groundPlane = 0.0;
	else
		groundPlane = 2.0 * sin(twoPiOverLambda * _dipoleSize * cosZenith);
	
	// normalize to zenith
	if(_zenithNorm)
		groundPlane /= 2.0 * sin(twoPiOverLambda * _dipoleSize);

	double
		sinDecAntennaZenith, cosDecAntennaZenith,
		sinDec, cosDec,
		sinHa, cosHa;
	sincos(decAntennaZenith, &sinDecAntennaZenith, &cosDecAntennaZenith);
	sincos(dec, &sinDec, &cosDec);
	sincos(ha - haAntennaZenith, &sinHa, &cosHa);
	
	// Notice that X and Y conventions are not equal to the RTS conventions.
	// Here, X is assumed to be East-West aligned (sensitive to EM radiation with
	// a East-West polarized vector).
	double rot[4];
	rot[0] =  cosHa;
	rot[1] =  sinDec*sinHa;
	rot[2] = -sinDecAntennaZenith*sinHa;
	rot[3] =  cosDecAntennaZenith*cosDec + sinDecAntennaZenith*sinDec*cosHa;
	//std::cout << "rot[0]=" << rot[0] << " groundPlane=" << groundPlane << " arrayFactor=" << arrayFactor << '\n';
	
	arrayFactor *= groundPlane;
	gain[0] = rot[0] * arrayFactor;
	gain[1] = rot[1] * arrayFactor;
	gain[2] = rot[2] * arrayFactor;
	gain[3] = rot[3] * arrayFactor;
}
