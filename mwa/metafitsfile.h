#ifndef METAFITSFILE_H
#define METAFITSFILE_H

#include <fitsio.h>

#include <vector>

#include "../fitsiochecker.h"

class MWAHeader;
class MWAHeaderExt;

class MetaFitsFile : private FitsIOChecker
{
public:
	MetaFitsFile(const char *filename);
	~MetaFitsFile();
	
	void ReadHeader(MWAHeader &header, MWAHeaderExt &headerExt);
	void ReadTiles(std::vector<class MWAInput> &inputs, std::vector<class MWAAntenna> &antennae);
	static void GetDelays(fitsfile* f, int* delays)
	{
		char keyValue[80];
		int status;
		fits_read_key(f, TSTRING, "DELAYS", keyValue, 0, &status);
		checkStatus(status, "GetDelays()");
		parseIntArray(keyValue, delays, 16);
	}
private:
	void parseKeyword(MWAHeader &header, MWAHeaderExt &headerExt, const std::string& keyName, const std::string& keyValue);
	static std::string parseFitsString(const char *valueStr);
	void parseFitsDate(const char *valueStr, int &year, int &month, int &day, int &hour, int &min, double &sec);
	static void parseIntArray(const char* valueStr, int* values, size_t count);
	bool parseBool(const char *valueStr);
	static std::string stripBand(const std::string &input);
	static bool isDigit(char c) { return c >= '0' && c <= '9'; }
	double parseFitsDateToMJD(const char* valueStr);
	
	std::string _filename;
	fitsfile *_fptr;
};

#endif
