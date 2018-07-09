#ifndef FITS_IO_CHECKER_H
#define FITS_IO_CHECKER_H

#include <string>

class FitsIOChecker
{
protected:
	static void checkStatus(int status, const std::string& filename);
	static void checkStatus(int status, const std::string& filename, const std::string& operation);
public:
	enum Unit {
		JanskyPerBeam,
		Jansky,
		Kelvin,
		MilliKelvin
	};
};

#endif
