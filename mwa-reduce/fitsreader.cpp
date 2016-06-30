#include "fitsreader.h"
#include "polarizationenum.h"

#include <stdexcept>
#include <sstream>
#include <cmath>

#include <casacore/fits/FITS/FITSDateUtil.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/measures/Measures/MeasConvert.h>

FitsReader::FitsReader(const FitsReader& source) :
	_filename(source._filename),
	_fitsPtr(0),
	_imgWidth(source._imgWidth), _imgHeight(source._imgHeight),
	_phaseCentreRA(source._phaseCentreRA), _phaseCentreDec(source._phaseCentreDec),
	_pixelSizeX(source._pixelSizeX), _pixelSizeY(source._pixelSizeY),
	_phaseCentreDL(source._phaseCentreDL), _phaseCentreDM(source._phaseCentreDM),
	_frequency(source._frequency), _bandwidth(source._bandwidth), _dateObs(source._dateObs),
	_hasBeam(source._hasBeam),
	_beamMajorAxisRad(source._beamMajorAxisRad), _beamMinorAxisRad(source._beamMinorAxisRad), _beamPositionAngle(source._beamPositionAngle),
	_polarization(source._polarization),
	_origin(source._origin), _originComment(source._originComment),
	_history(source._history)
{
	int status = 0;
	fits_open_file(&_fitsPtr, _filename.c_str(), READONLY, &status);
	checkStatus(status, _filename);
	
	// Move to first HDU
	int hduType;
	fits_movabs_hdu(_fitsPtr, 1, &hduType, &status);
	checkStatus(status, _filename);
	if(hduType != IMAGE_HDU) throw std::runtime_error("First HDU is not an image");
}

FitsReader::~FitsReader()
{
	int status = 0;
	fits_close_file(_fitsPtr, &status);
}

FitsReader& FitsReader::operator=(const FitsReader& rhs)
{
	_filename = rhs._filename;
	_imgWidth = rhs._imgWidth;
	_imgHeight = rhs._imgHeight;
	_phaseCentreRA = rhs._phaseCentreRA;
	_phaseCentreDec = rhs._phaseCentreDec;
	_pixelSizeX = rhs._pixelSizeX;
	_pixelSizeY = rhs._pixelSizeY;
	_phaseCentreDL = rhs._phaseCentreDL;
	_phaseCentreDM = rhs._phaseCentreDM;
	_frequency = rhs._frequency;
	_bandwidth = rhs._bandwidth;
	_dateObs = rhs._dateObs;
	_hasBeam = rhs._hasBeam;
	_beamMajorAxisRad = rhs._beamMajorAxisRad;
	_beamMinorAxisRad = rhs._beamMinorAxisRad;
	_beamPositionAngle = rhs._beamPositionAngle;
	_polarization = rhs._polarization;
	_origin = rhs._origin;
	_originComment = rhs._originComment;
	_history = rhs._history;
	
	int status = 0;
	fits_close_file(_fitsPtr, &status);
	checkStatus(status, _filename);
	
	fits_open_file(&_fitsPtr, _filename.c_str(), READONLY, &status);
	checkStatus(status, _filename);
	
	// Move to first HDU
	int hduType;
	fits_movabs_hdu(_fitsPtr, 1, &hduType, &status);
	checkStatus(status, _filename);
	if(hduType != IMAGE_HDU) throw std::runtime_error("First HDU is not an image");
	return *this;
}

float FitsReader::readFloatKey(const char *key)
{
	int status = 0;
	float floatValue;
	fits_read_key(_fitsPtr, TFLOAT, key, &floatValue, 0, &status);
	checkStatus(status, _filename, std::string("Read float key ") + key);
	return floatValue;
}

double FitsReader::readDoubleKey(const char *key)
{
	int status = 0;
	double value;
	fits_read_key(_fitsPtr, TDOUBLE, key, &value, 0, &status);
	checkStatus(status, _filename, std::string("Read float key ") + key);
	return value;
}

bool FitsReader::ReadFloatKeyIfExists(const char *key, float &dest)
{
	int status = 0;
	float floatValue;
	fits_read_key(_fitsPtr, TFLOAT, key, &floatValue, 0, &status);
	if(status == 0)
		dest = floatValue;
	return status == 0;
}

bool FitsReader::ReadDoubleKeyIfExists(const char *key, double &dest)
{
	int status = 0;
	double doubleValue;
	fits_read_key(_fitsPtr, TDOUBLE, key, &doubleValue, 0, &status);
	if(status == 0)
		dest = doubleValue;
	return status == 0;
}

bool FitsReader::readDateKeyIfExists(const char *key, double &dest)
{
	int status = 0;
	char keyStr[256];
	fits_read_key(_fitsPtr, TSTRING, key, keyStr, 0, &status);
	if(status == 0)
	{
		dest = FitsReader::ParseFitsDateToMJD(keyStr);
		return true;
	}
	else return false;
}

std::string FitsReader::readStringKey(const char *key)
{
	int status = 0;
	char keyStr[256];
	fits_read_key(_fitsPtr, TSTRING, key, keyStr, 0, &status);
	checkStatus(status, _filename, std::string("Read string key ") + key);
	return std::string(keyStr);
}

bool FitsReader::ReadStringKeyIfExists(const char *key, std::string& value, std::string& comment)
{
	int status = 0;
	char valueStr[256], commentStr[256];
	fits_read_key(_fitsPtr, TSTRING, key, valueStr, commentStr, &status);
	if(status == 0)
	{
		value = valueStr;
		comment = commentStr;
	}
	return status == 0;
}

void FitsReader::initialize()
{
	int status = 0;
	fits_open_file(&_fitsPtr, _filename.c_str(), READONLY, &status);
	checkStatus(status, _filename);
	
	// Move to first HDU
	int hduType;
	fits_movabs_hdu(_fitsPtr, 1, &hduType, &status);
	checkStatus(status, _filename);
	if(hduType != IMAGE_HDU) throw std::runtime_error("First HDU is not an image");
	
	int naxis = 0;
	fits_get_img_dim(_fitsPtr, &naxis, &status);
	checkStatus(status, _filename);
	if(naxis < 2) throw std::runtime_error("NAxis in image < 2");
	
	std::vector<long> naxes(naxis);
	fits_get_img_size(_fitsPtr, naxis, &naxes[0], &status);
	checkStatus(status, _filename);
	for(int i=2;i!=naxis;++i)
		if(naxes[i] != 1) throw std::runtime_error("Multiple images in fits file");
	_imgWidth = naxes[0];
	_imgHeight = naxes[1];
	
	double bScale = 1.0, bZero = 0.0, equinox = 2000.0;
	ReadDoubleKeyIfExists("BSCALE", bScale);
	ReadDoubleKeyIfExists("BZERO", bZero);
	ReadDoubleKeyIfExists("EQUINOX", equinox);
	if(bScale != 1.0)
		throw std::runtime_error("Invalid value for BSCALE");
	if(bZero != 0.0)
		throw std::runtime_error("Invalid value for BZERO");
	if(equinox != 2000.0)
		throw std::runtime_error("Invalid value for EQUINOX: "+readStringKey("EQUINOX"));
	
	std::string tmp;
	if(ReadStringKeyIfExists("CTYPE1", tmp) && tmp != "RA---SIN")
		throw std::runtime_error("Invalid value for CTYPE1");
	
	_phaseCentreRA = 0.0;
	ReadDoubleKeyIfExists("CRVAL1", _phaseCentreRA);
	_phaseCentreRA *= M_PI / 180.0;
	_pixelSizeX = 0.0;
	ReadDoubleKeyIfExists("CDELT1", _pixelSizeX);
	_pixelSizeX *= -M_PI / 180.0;
	if(ReadStringKeyIfExists("CUNIT1", tmp) && tmp != "deg")
		throw std::runtime_error("Invalid value for CUNIT1");
	double centrePixelX = 0.0;
	if(ReadDoubleKeyIfExists("CRPIX1", centrePixelX))
		_phaseCentreDL = (centrePixelX - ((_imgWidth / 2.0)+1.0)) * _pixelSizeX;
	else
		_phaseCentreDL = 0.0;

	if(ReadStringKeyIfExists("CTYPE2",tmp) && tmp != "DEC--SIN")
		throw std::runtime_error("Invalid value for CTYPE2");
	_phaseCentreDec = 0.0;
	ReadDoubleKeyIfExists("CRVAL2", _phaseCentreDec);
	_phaseCentreDec *= M_PI / 180.0;
	_pixelSizeY = 0.0;
	ReadDoubleKeyIfExists("CDELT2", _pixelSizeY);
	_pixelSizeY *= M_PI / 180.0;
	if(ReadStringKeyIfExists("CUNIT2", tmp) && tmp != "deg")
		throw std::runtime_error("Invalid value for CUNIT2");
	double centrePixelY = 0.0;
	if(ReadDoubleKeyIfExists("CRPIX2", centrePixelY))
		_phaseCentreDM = ((_imgHeight / 2.0)+1.0 - centrePixelY) * _pixelSizeY;
	else
		_phaseCentreDM = 0.0;
	
	_dateObs = 0.0;
	readDateKeyIfExists("DATE-OBS", _dateObs);
	
	if(naxis >= 3 && readStringKey("CTYPE3") == "FREQ")
	{
		_frequency = readDoubleKey("CRVAL3");
		_bandwidth = readDoubleKey("CDELT3");
	}
	else {
		_frequency = 0.0;
		_bandwidth = 0.0;
	}
	
	if(naxis >= 4 && readStringKey("CTYPE4") == "STOKES")
	{
		double val = readDoubleKey("CRVAL4");
		switch(int(val))
		{
			default: throw std::runtime_error("Unknown polarization specified in fits file");
			case 1: _polarization = Polarization::StokesI; break;
			case 2: _polarization = Polarization::StokesQ; break;
			case 3: _polarization = Polarization::StokesU; break;
			case 4: _polarization = Polarization::StokesV; break;
			case -1: _polarization = Polarization::RR; break;
			case -2: _polarization = Polarization::LL; break;
			case -3: _polarization = Polarization::RL; break;
			case -4: _polarization = Polarization::LR; break;
			case -5: _polarization = Polarization::XX; break;
			case -6: _polarization = Polarization::YY; break;
			case -7: _polarization = Polarization::XY; break;
			case -8: _polarization = Polarization::YX; break;
		}
	}
	else {
		_polarization = Polarization::StokesI;
	}
	double bMaj=0.0, bMin=0.0, bPa=0.0;
	if(ReadDoubleKeyIfExists("BMAJ", bMaj) && ReadDoubleKeyIfExists("BMIN", bMin) && ReadDoubleKeyIfExists("BPA", bPa))
	{
		_hasBeam = true;
		_beamMajorAxisRad = bMaj * (M_PI / 180.0);
		_beamMinorAxisRad = bMin * (M_PI / 180.0);
		_beamPositionAngle = bPa * (M_PI / 180.0);
	}
	else {
		_hasBeam = false;
		_beamMajorAxisRad = 0.0;
		_beamMinorAxisRad = 0.0;
		_beamPositionAngle = 0.0;
	}
	
	_origin = std::string();
	_originComment = std::string();
	ReadStringKeyIfExists("ORIGIN", _origin, _originComment);
	
	_history.clear();
	readHistory();
}

template void FitsReader::Read(float* image);
template void FitsReader::Read(double* image);

template<typename NumType>
void FitsReader::Read(NumType* image)
{
	int status = 0;
	int naxis = 0;
	fits_get_img_dim(_fitsPtr, &naxis, &status);
	checkStatus(status, _filename);
	std::vector<long> firstPixel(naxis);
	for(int i=0;i!=naxis;++i) firstPixel[i] = 1;
	
	if(sizeof(NumType)==8)
		fits_read_pix(_fitsPtr, TDOUBLE, &firstPixel[0], _imgWidth*_imgHeight, 0, image, 0, &status);
	else if(sizeof(NumType)==4)
		fits_read_pix(_fitsPtr, TFLOAT, &firstPixel[0], _imgWidth*_imgHeight, 0, image, 0, &status);
	else
		throw std::runtime_error("sizeof(NumType)!=8 || 4 not implemented");
	checkStatus(status, _filename);
}

void FitsReader::readHistory()
{
	int status = 0;
	int npos, moreKeys;
	fits_get_hdrspace(_fitsPtr, &npos, &moreKeys, &status);
	checkStatus(status, _filename);
	char keyCard[256];
	for(int pos=1; pos<=npos; ++pos)
	{
		fits_read_record(_fitsPtr, pos, keyCard, &status);
		keyCard[7] = 0;
		if(std::string("HISTORY") == keyCard) {
			_history.push_back(&keyCard[8]);
		}
	}
}

double FitsReader::ParseFitsDateToMJD(const char* valueStr)
{
	casacore::MVTime time;
	casacore::MEpoch::Types systypes;
	bool parseSuccess = casacore::FITSDateUtil::fromFITS(time, systypes, valueStr, "UTC");
	if(!parseSuccess)
		throw std::runtime_error(std::string("Could not parse FITS date: ") + valueStr);
	casacore::MEpoch epoch(time.get(), systypes);
	return epoch.getValue().get();
}
