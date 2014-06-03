#ifndef JONES_LOOKUP_DIPOLE_H
#define JONES_LOOKUP_DIPOLE_H

#include "../fitsiochecker.h"

#include "../aocommon/uvector.h"

#include <complex>
#include <map>

#include "../alglib/ap.h"
#include "../alglib/interpolation.h"

/**
 * This class is based on Python code written by Randall Wayth,
 * function getJonesLookup(self,az,za,freq,zenith_norm=True) in
 * mwapy/pb/mwa_tile.py, 2014-06-02.
 */
class JonesLookupDipole : private FitsIOChecker
{
public:
	static void Initialize() { if(_tables.empty()) loadLookupTable(); }
	
	static void Interpolate(std::complex<double>* jonesMatrix, double az, double za, double freq, bool zenithNorm = true);
private:
	struct FrequencyTable
	{
		double frequency;
		alglib::real_1d_array values;
		alglib::spline2dinterpolant spline;
		std::complex<double> norm[4];
	};
	
	static void loadLookupTable();
	
	static std::complex<double>& d2c(double& v)
	{
		return *reinterpret_cast<std::complex<double>*>(&v);
	}
	
	static alglib::real_1d_array _zaValues, _phValues;
	
	static std::map<double, FrequencyTable> _tables;
};

#endif // JONES_LOOKUP_DIPOLE_H
