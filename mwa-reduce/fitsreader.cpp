#include "fitsreader.h"

#include <stdexcept>
#include <sstream>
#include <cmath>

void FitsReader::checkStatus(int status, const std::string &operation) 
{
	if(status) {
		/* fits_get_errstatus returns at most 30 characters */
		char err_text[31];
		fits_get_errstatus(status, err_text);
		char err_msg[81];
		std::stringstream errMsg;
		if(!operation.empty())
			errMsg << "During operation " << operation << ", ";
		errMsg << "CFITSIO reported error when performing IO on file '" << _filename << "': " << err_text << " (";
		while(fits_read_errmsg(err_msg))
			errMsg << err_msg;
		errMsg << ')';
		throw std::runtime_error(errMsg.str());
	}
}

float FitsReader::readFloatKey(const char *key)
{
	int status = 0;
	float floatValue;
	fits_read_key(_fitsPtr, TFLOAT, key, &floatValue, 0, &status);
	checkStatus(status, key);
	return floatValue;
}

void FitsReader::readFloatKeyIfExists(const char *key, float &dest)
{
	int status = 0;
	float floatValue;
	fits_read_key(_fitsPtr, TFLOAT, key, &floatValue, 0, &status);
	if(status == 0)
		dest = floatValue;
}

std::string FitsReader::readStringKey(const char *key)
{
	int status = 0;
	char keyStr[256];
	fits_read_key(_fitsPtr, TSTRING, key, keyStr, 0, &status);
	checkStatus(status, key);
	return std::string(keyStr);
}

void FitsReader::initialize()
{
	int status = 0;
	fits_open_file(&_fitsPtr, _filename.c_str(), READONLY, &status);
	checkStatus(status);
	
	// Move to first HDU
	int hduType;
	fits_movabs_hdu(_fitsPtr, 1, &hduType, &status);
	checkStatus(status);
	if(hduType != IMAGE_HDU) throw std::runtime_error("First HDU is not an image");
	
	int naxis = 0;
	fits_get_img_dim(_fitsPtr, &naxis, &status);
	checkStatus(status);
	if(naxis < 2) throw std::runtime_error("NAxis in image < 2");
	
	long naxes[naxis];
	fits_get_img_size(_fitsPtr, naxis, naxes, &status);
	checkStatus(status);
	for(int i=2;i!=naxis;++i)
		if(naxes[i] != 1) throw std::runtime_error("Multiple images in fits file");
	_imgWidth = naxes[0];
	_imgHeight = naxes[1];
	
	float bScale = 1.0, bZero = 0.0, equinox = 2000.0;
	readFloatKeyIfExists("BSCALE", bScale);
	readFloatKeyIfExists("BZERO", bZero);
	readFloatKeyIfExists("EQUINOX", equinox);
	if(bScale != 1.0)
		throw std::runtime_error("Invalid value for BSCALE");
	if(bZero != 0.0)
		throw std::runtime_error("Invalid value for BZERO");
	if(equinox != 2000.0)
		throw std::runtime_error("Invalid value for EQUINOX: "+readStringKey("EQUINOX"));
	
	if(readStringKey("CTYPE1") != "RA---SIN")
		throw std::runtime_error("Invalid value for CTYPE1");
	_phaseCentreRA = readFloatKey("CRVAL1") * (M_PI / 180.0);
	_pixelSizeX = readFloatKey("CDELT1") * (-M_PI / 180.0);
	if(readStringKey("CUNIT1") != "deg")
		throw std::runtime_error("Invalid value for CUNIT1");
	
	if(readStringKey("CTYPE2") != "DEC--SIN")
		throw std::runtime_error("Invalid value for CTYPE2");
	_phaseCentreDec = readFloatKey("CRVAL2") * (M_PI / 180.0);
	_pixelSizeY = readFloatKey("CDELT2") * (M_PI / 180.0);
	if(readStringKey("CUNIT2") != "deg")
		throw std::runtime_error("Invalid value for CUNIT2");
}

template void FitsReader::Read(double* image);

template<typename NumType>
void FitsReader::Read(NumType* image)
{
	int status = 0;
	int naxis = 0;
	fits_get_img_dim(_fitsPtr, &naxis, &status);
	checkStatus(status);
	long firstPixel[naxis];
	for(int i=0;i!=naxis;++i) firstPixel[i] = 1;
	
	if(sizeof(NumType)==8)
		fits_read_pix(_fitsPtr, TDOUBLE, firstPixel, _imgWidth*_imgHeight, 0, image, 0, &status);
	else
		throw std::runtime_error("sizeof(NumType)!=8 not implemented");
	checkStatus(status);
}

FitsReader::~FitsReader()
{
	int status = 0;
	fits_close_file(_fitsPtr, &status);
	checkStatus(status);
}
