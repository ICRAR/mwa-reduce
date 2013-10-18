#ifndef FITSWRITER_H
#define FITSWRITER_H

#include <string>

#include "polarizationenum.h"

class FitsWriter
{
	public:
		FitsWriter() :
			_width(0), _height(0),
			_phaseCentreRA(0.0), _phaseCentreDec(0.0), _pixelSizeX(0.0), _pixelSizeY(0.0),
			_frequency(0.0), _bandwidth(0.0),
			_dateObs(0.0),
			_hasBeam(false),
			_beamMajorAxisRad(0.0), _beamMinorAxisRad(0.0), _beamPositionAngle(0.0),
			_polarization(Polarization::StokesI)
		{
		}
		
		FitsWriter(const class FitsReader& reader) :
			_width(0), _height(0),
			_phaseCentreRA(0.0), _phaseCentreDec(0.0), _pixelSizeX(0.0), _pixelSizeY(0.0),
			_frequency(0.0), _bandwidth(0.0),
			_dateObs(0.0),
			_hasBeam(false),
			_beamMajorAxisRad(0.0), _beamMinorAxisRad(0.0), _beamPositionAngle(0.0),
			_polarization(Polarization::StokesI)
		{
			SetMetadata(reader);
		}
		
		template<typename NumType> void Write(const std::string& filename, const NumType* image);
		
		void SetBeamInfo(double widthRad)
		{
			SetBeamInfo(widthRad, widthRad, 0.0);
		}
		void SetBeamInfo(double majorAxisRad, double minorAxisRad, double positionAngleRad)
		{
			_hasBeam = true;
			_beamMajorAxisRad = majorAxisRad;
			_beamMinorAxisRad = minorAxisRad;
			_beamPositionAngle = positionAngleRad;
		}
		void SetImageDimensions(size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY)
		{
			_width = width;
			_height = height;
			_phaseCentreRA = phaseCentreRA;
			_phaseCentreDec = phaseCentreDec;
			_pixelSizeX = pixelSizeX;
			_pixelSizeY = pixelSizeY;
		}
		void SetFrequency(double frequency, double bandwidth)
		{
			_frequency = frequency;
			_bandwidth = bandwidth;
		}
		void SetDate(double dateObs)
		{
			_dateObs = dateObs;
		}
		void SetPolarization(PolarizationEnum polarization)
		{
			_polarization = polarization;
		}
		void SetMetadata(const class FitsReader& reader);
		
		double RA() const { return _phaseCentreRA; }
		double Dec() const { return _phaseCentreDec; }
		double Frequency() const { return _frequency; }
		double Bandwidth() const { return _bandwidth; }
		double BeamSizeMajorAxis() const { return _beamMajorAxisRad; }
	private:
		std::size_t _width, _height;
		double _phaseCentreRA, _phaseCentreDec, _pixelSizeX, _pixelSizeY;
		double _frequency, _bandwidth;
		double _dateObs;
		bool _hasBeam;
		double _beamMajorAxisRad, _beamMinorAxisRad, _beamPositionAngle;
		PolarizationEnum _polarization;
		
		void checkStatus(int status, const std::string& filename);
		void julianDateToYMD(double jd, int &year, int &month, int &day);
		void mjdToHMS(double mjd, int& hour, int& minutes, int& seconds, int& deciSec);
};

#endif
