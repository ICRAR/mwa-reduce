#include "fitswriter.h"

#include <stdexcept>
#include <sstream>
#include <vector>

#include <fitsio2.h>
#include <cmath>

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
void FitsWriter::Write(NumType *image, size_t imgSize, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY)
{
	int status = 0;
	fitsfile *fptr;
	fits_create_file(&fptr, (std::string("!") + _filename).c_str(), &status);
	checkStatus(status);
	
	// append image HDU
	int bitPixInt = FLOAT_IMG;
	long naxes[2];
	naxes[0] = imgSize;
	naxes[1] = imgSize;
	fits_create_img(fptr, bitPixInt, 2, naxes, &status);
	checkStatus(status);
	float zero = 0, one = 1, equinox = 2000.0;
	fits_write_key(fptr, TFLOAT, "BSCALE", (void*) &one, "", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "BZERO", (void*) &zero, "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "BUNIT", (void*) "JY/BEAM", "Units are in Jansky per beam", &status); checkStatus(status);
	fits_write_key(fptr, TFLOAT, "EQUINOX", (void*) &equinox, "J2000", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "BTYPE", (void*) "Intensity", "", &status); checkStatus(status);
	fits_write_key(fptr, TSTRING, "ORIGIN", (void*) "AOImager", "Imager written by Andre Offringa", &status); checkStatus(status);
	float phaseCentreRADeg = (phaseCentreRA/M_PI)*180.0, phaseCentreDecDeg = (phaseCentreDec/M_PI)*180.0;
	float centrePixelX = (imgSize / 2.0)+1, centrePixelY = (imgSize / 2.0)+1;
	float stepXDeg = -(pixelSizeX/M_PI)*180.0, stepYDeg = (pixelSizeY/M_PI)*180.0;
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
	// RESTFRQ ?
	fits_write_key(fptr, TSTRING, "SPECSYS", (void*) "TOPOCENT", "", &status); checkStatus(status);
	
	long firstpixel[2];
	for(int i=0;i < 2;i++) firstpixel[i] = 1;
	if(sizeof(NumType)==8)
	{
		double nullValue = 0.0;
		fits_write_pixnull(fptr, TDOUBLE, firstpixel, imgSize*imgSize, image, &nullValue, &status);
	}
	else if(sizeof(NumType)==4)
	{
		float nullValue = 0.0;
		fits_write_pixnull(fptr, TFLOAT, firstpixel, imgSize*imgSize, image, &nullValue, &status);
	}
	else
	{
		double nullValue = 0.0;
		size_t totalSize = imgSize * imgSize;
		std::vector<double> copy(totalSize);
		for(size_t i=0;i!=totalSize;++i) copy[i] = image[i];
		fits_write_pixnull(fptr, TDOUBLE, firstpixel, imgSize*imgSize, &copy[0], &nullValue, &status);
	}
	checkStatus(status);
	
	fits_close_file(fptr, &status);
	checkStatus(status);
}

template void FitsWriter::Write<long double>(long double *image, size_t imgSize, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY);
template void FitsWriter::Write<double>(double *image, size_t imgSize, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY);
template void FitsWriter::Write<float>(float *image, size_t imgSize, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY);
