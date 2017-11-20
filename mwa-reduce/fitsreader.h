#ifndef FITSREADER_H
#define FITSREADER_H

#include <string>
#include <vector>

#include <fitsio.h>

#include "polarization.h"
#include "fitsiochecker.h"

class FitsReader : public FitsIOChecker
{
	public:
		explicit FitsReader(const std::string &filename) 
		: _filename(filename), _hasBeam(false)
		{
			initialize(); 
		}
		FitsReader(const FitsReader& source);
		~FitsReader();
		
		FitsReader& operator=(const FitsReader& rhs);
		
		template<typename NumType> void Read(NumType *image);
		
		size_t ImageWidth() const { return _imgWidth; }
		size_t ImageHeight() const { return _imgHeight; }
		
		double PhaseCentreRA() const { return _phaseCentreRA; }
		double PhaseCentreDec() const { return _phaseCentreDec; }
		
		double PixelSizeX() const { return _pixelSizeX; }
		double PixelSizeY() const { return _pixelSizeY; }
		
		double PhaseCentreDL() const { return _phaseCentreDL; }
		double PhaseCentreDM() const { return _phaseCentreDM; }
		
		double Frequency() const { return _frequency; }
		double Bandwidth() const { return _bandwidth; }
		
		double DateObs() const { return _dateObs; }
		PolarizationEnum Polarization() const { return _polarization; }
		
		FitsIOChecker::Unit Unit() const { return _unit; }
		
		bool HasBeam() const { return _hasBeam; }
		double BeamMajorAxisRad() const { return _beamMajorAxisRad; }
		double BeamMinorAxisRad() const { return _beamMinorAxisRad; }
		double BeamPositionAngle() const { return _beamPositionAngle; }
		
		const std::string& TelescopeName() const { return _telescopeName; }
		const std::string& Observer() const { return _observer; }
		const std::string& ObjectName() const { return _objectName; }
		
		const std::string& Origin() const { return _origin; }
		const std::string& OriginComment() const { return _originComment; }
		
		const std::vector<std::string>& History() const { return _history; }
		
		bool ReadDoubleKeyIfExists(const char* key, double& dest);
		bool ReadStringKeyIfExists(const char* key, std::string& dest) {
			std::string c;
			return ReadStringKeyIfExists(key, dest, c);
		}
		bool ReadStringKeyIfExists(const char* key, std::string& value, std::string& comment);
		bool ReadFloatKeyIfExists(const char* key, float& dest);
		
		static double ParseFitsDateToMJD(const char* valueStr);
		
		const std::string& Filename() const { return _filename; }
		
		fitsfile* FitsHandle() const { return _fitsPtr; }
	private:
		double readDoubleKey(const char* key);
		std::string readStringKey(const char* key);
		void readHistory();
		bool readDateKeyIfExists(const char *key, double &dest);
		
		void initialize();
		
		std::string _filename;
		fitsfile *_fitsPtr;
		
		size_t _imgWidth, _imgHeight;
		double _phaseCentreRA, _phaseCentreDec;
		double _pixelSizeX, _pixelSizeY;
		double _phaseCentreDL, _phaseCentreDM;
		double _frequency, _bandwidth, _dateObs;
		bool _hasBeam;
		double _beamMajorAxisRad, _beamMinorAxisRad, _beamPositionAngle;
		
		PolarizationEnum _polarization;
		FitsIOChecker::Unit _unit;
		std::string _telescopeName, _observer, _objectName;
		std::string _origin, _originComment;
		std::vector<std::string> _history;
};

#endif
