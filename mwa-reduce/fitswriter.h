#ifndef FITSWRITER_H
#define FITSWRITER_H

#include <string>

class FitsWriter
{
	public:
		FitsWriter(const std::string &filename) 
		: _filename(filename) { }
		template<typename NumType> void Write(const NumType *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs);
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
	private:
		std::string _filename;
		bool _hasBeam;
		double _beamMajorAxisRad, _beamMinorAxisRad, _beamPositionAngle;
		
		void checkStatus(int status);
		void julianDateToYMD(double jd, int &year, int &month, int &day);
};

#endif
