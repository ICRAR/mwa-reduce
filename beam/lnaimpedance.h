#ifndef LNA_IMPEDANCE_H
#define LNA_IMPEDANCE_H

#include <cstring>
#include <complex>

class LNAImpedance
{
	public:
		static std::complex<double> Get(double frequency)
		{
			double mhz = frequency*1e-6;
			size_t indexLow = size_t(mhz)-50;
			// Do linear interpolation
			double modulo = mhz - floor(mhz);
			return impedanceArray[indexLow]*(1.0-modulo) + impedanceArray[indexLow+1]*modulo;
		}
	private:
		typedef std::complex<double> ctype;
		static double impedanceFrequency(size_t index)
		{
			return (50.0+double(index)*1e6);
		}
		
    /** Measured MWA LNA impedance between 50 and 500 MHz.
		 */
		static const ctype impedanceArray[452];
};

#endif
