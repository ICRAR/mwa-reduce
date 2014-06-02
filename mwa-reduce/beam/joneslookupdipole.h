#ifndef JONES_LOOKUP_DIPOLE_H
#define JONES_LOOKUP_DIPOLE_H

#ifdef HAVE_ALGLIB

#include "../fitsiochecker.h"

#include "../aocommon/uvector.h"

#include <complex>
#include <map>

#include <ap.h>
#include <interpolation.h>

/**
 * This class is based on Python code written by Randall Wayth,
 * function getJonesLookup(self,az,za,freq,zenith_norm=True) in
 * mwapy/pb/mwa_tile.py, 2014-06-02.
 */
class JonesLookupDipole : private FitsIOChecker
{
public:
	static void Initialize() { if(_tables.empty()) loadLookupTable("Jmatrix.fits"); }
	
	static void Interpolate(std::complex<double>* jonesMatrix, double az, double za, double freq, bool zenithNorm = true);
private:
	struct FrequencyTable
	{
		double frequency;
		alglib::real_1d_array values;
		alglib::spline2dinterpolant spline;
		std::complex<double> norm[4];
	};
	
	static void loadLookupTable(const std::string& filename);
	
	static std::complex<double>& d2c(double& v)
	{
		return *reinterpret_cast<std::complex<double>*>(&v);
	}
	
	static alglib::real_1d_array _zaValues, _phValues;
	
	static std::map<double, FrequencyTable> _tables;
};

#else // HAVE_ALGLIB

#include <complex>
#include <stdexcept>

class JonesLookupDipole
{
public:
	static void Initialize()
	{
		throw std::runtime_error("Can not evaluate Jones look-up dipole: library alglib missing");
	}
	
	static void Interpolate(std::complex<double>* jonesMatrix, double az, double za, double freq, bool zenithNorm = true)
	{
		throw std::runtime_error("Can not evaluate Jones look-up dipole: library alglib missing");
	}
};

#endif // HAVE_ALGLIB

#endif // JONES_LOOKUP_DIPOLE_H
