#include "fitswriter.h"

#include <stdexcept>
#include <sstream>
#include <vector>

#include <fitsio2.h>
#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>

void FitsWriter::checkStatus(int status) 
{
	if(status) {
		/* fits_get_errstatus returns at most 30 characters */
		char err_text[31];
		fits_get_errstatus(status, err_text);
		char err_msg[81];
		std::stringstream errMsg;
		errMsg << "CFITSIO reported error when performing IO on file '" << _filename << "':" << err_text << " (";
		while(fits_read_errmsg(err_msg))
			errMsg << err_msg;
		errMsg << ')';
		throw std::runtime_error(errMsg.str());
	}
}

template<typename NumType>
void FitsWriter::Write(const NumType *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs)
{
	int status = 0;
	fitsfile *fptr;
	fits_create_file(&fptr, (std::string("!") + _filename).c_str(), &status);
	checkStatus(status);
	
	// append image HDU
	int bitPixInt = FLOAT_IMG;
	long naxes[4];
	naxes[0] = width;
	naxes[1] = height;
	naxes[2] = 1;
	naxes[3] = 1;
	fits_create_img(fptr, bitPixInt, 4, naxes, &status);
	checkStatus(status);
	float zero = 0, one = 1, equinox = 2000.0;
	fits_write_key(fptr, TFLOAT, "BSCALE", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "BZERO", (void*) &zero, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "BUNIT", (void*) "JY/BEAM", "Units are in Jansky per beam", &status); checkStatus(status);
	
	if(_hasBeam)
	{
		float
			majDeg = _beamMajorAxisRad * 180.0 / M_PI,
			minDeg = _beamMinorAxisRad * 180.0 / M_PI, 
			posAngle = _beamPositionAngle * 180.0 / M_PI;
		fits_write_key(fptr, TFLOAT, "BMAJ", (void*) &majDeg, "", &status); checkStatus(status);
		fits_write_key(fptr, TFLOAT, "BMIN", (void*) &minDeg, "", &status); checkStatus(status);
		fits_write_key(fptr, TFLOAT, "BPA", (void*) &posAngle, "", &status); checkStatus(status);
	}
	
	fits_write_key(fptr, TFLOAT, "EQUINOX", (void*) &equinox, "J2000", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "BTYPE", (void*) "Intensity", "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "ORIGIN", (void*) "AO/WSImager", "Imager written by Andre Offringa", &status); checkStatus(status);
	float phaseCentreRADeg = (phaseCentreRA/M_PI)*180.0, phaseCentreDecDeg = (phaseCentreDec/M_PI)*180.0;
	float centrePixelX = (width / 2.0)+1, centrePixelY = (height / 2.0)+1;
	float stepXDeg = (-pixelSizeX/M_PI)*180.0, stepYDeg = (pixelSizeY/M_PI)*180.0;
	fits_write_key(fptr, TSTRING, "CTYPE1", (void*) "RA---SIN", "Right ascension angle cosine", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRPIX1", (void*) &centrePixelX, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRVAL1", (void*) &phaseCentreRADeg, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CDELT1", (void*) &stepXDeg, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "CUNIT1", (void*) "deg", "", &status); checkStatus(status);
	
	fits_write_key(fptr, TSTRING, "CTYPE2", (void*) "DEC--SIN", "Declination angle cosine", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRPIX2", (void*) &centrePixelY, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRVAL2", (void*) &phaseCentreDecDeg, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CDELT2", (void*) &stepYDeg, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "CUNIT2", (void*) "deg", "", &status); checkStatus(status);

	fits_write_key(fptr, TSTRING, "CTYPE3", (void*) "FREQ", "Central frequency", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRPIX3", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TDOUBLE, "CRVAL3", (void*) &frequency, "", &status); checkStatus(status);
	fits_write_key(fptr, TDOUBLE, "CDELT3", (void*) &bandwidth, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "CUNIT3", (void*) "Hz", "", &status); checkStatus(status);
	
	fits_write_key(fptr, TSTRING, "CTYPE4", (void*) "STOKES", "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRPIX4", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CRVAL4", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "CDELT4", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "CUNIT4", (void*) "Hz", "", &status); checkStatus(status);
	
	// RESTFRQ ?
	fits_write_key(fptr, TSTRING, "SPECSYS", (void*) "TOPOCENT", "", &status); checkStatus(status);
	
  int year, month, day, hour, min, sec, deciSec;
	std::cout << "Date obs= " << dateObs << '\n';
	julianDateToYMD(dateObs + 2400000.5, year, month, day);
	mjdToHMS(dateObs, hour, min, sec, deciSec);
	char dateStr[40];
  std::sprintf(dateStr, "%d-%02d-%02dT%02d:%02d:%02d.%01d", year, month, day, hour, min, sec, deciSec);
	fits_write_key(fptr, TSTRING, "DATE-OBS", (void*) dateStr, "", &status); checkStatus(status);
	
	long firstpixel[4];
	for(int i=0;i < 4;i++) firstpixel[i] = 1;
	if(sizeof(NumType)==8)
	{
		double nullValue = std::numeric_limits<double>::max();
		fits_write_pixnull(fptr, TDOUBLE, firstpixel, width*height, const_cast<double*>(reinterpret_cast<const double*>(image)), &nullValue, &status);
	}
	else if(sizeof(NumType)==4)
	{
		float nullValue = std::numeric_limits<float>::max();
		fits_write_pixnull(fptr, TFLOAT, firstpixel, width*height, const_cast<float*>(reinterpret_cast<const float*>(image)), &nullValue, &status);
	}
	else
	{
		double nullValue = std::numeric_limits<double>::max();
		size_t totalSize = width*height;
		std::vector<double> copy(totalSize);
		for(size_t i=0;i!=totalSize;++i) copy[i] = image[i];
		fits_write_pixnull(fptr, TDOUBLE, firstpixel, totalSize, &copy[0], &nullValue, &status);
	}
	checkStatus(status);
	
	fits_close_file(fptr, &status);
	checkStatus(status);
}

template void FitsWriter::Write<long double>(const long double *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs);
template void FitsWriter::Write<double>(const double *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs);
template void FitsWriter::Write<float>(const float *image, size_t width, size_t height, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double frequency, double bandwidth, double dateObs);

void FitsWriter::julianDateToYMD(double jd, int &year, int &month, int &day)
{
  int z = jd+0.5;
  int w = (z-1867216.25)/36524.25;
  int x = w/4;
  int a = z+1+w-x;
  int b = a+1524;
  int c = (b-122.1)/365.25;
  int d = 365.25*c;
  int e = (b-d)/30.6001;
  int f = 30.6001*e;
  day = b-d-f;
  while (e-1 > 12) e-=12;
  month = e-1;
  year = c-4715-((e-1)>2?1:0);
}

void FitsWriter::mjdToHMS(double mjd, int& hour, int& minutes, int& seconds, int& deciSec)
{
	hour = int(fmod(mjd * 24.0, 24.0));
	minutes = int(fmod(mjd*60.0 * 24.0, 60.0));
	seconds = int(fmod(mjd*3600.0 * 24.0, 60.0));
	deciSec = int(fmod(mjd*36000.0 * 24.0, 10.0));
}
