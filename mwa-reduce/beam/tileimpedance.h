#ifndef TILE_IMPEDANCE_H
#define TILE_IMPEDANCE_H

#include <complex>
#include <iostream>
#include <cstring>
#include <vector>

#define TILE_IMPEDANCE_MATRIX_COUNT 5

class TileImpedance
{
public:
	static const std::complex<double> *Get(double frequency)
	{
		if(matrices.empty())
			initializeMatrices();
		
		const struct ImpedanceMatrix *m = 0;
		double minDist = 1e10;
		
		for(size_t i=0; i!=TILE_IMPEDANCE_MATRIX_COUNT; ++i)
		{
			double d = fabs(matrices[i]._frequency - frequency);
			if(d < minDist) {
				minDist = d;
				m = &matrices[i];
			}
		}
		if(minDist > 2e6)
			std::cerr << "Nearest tabulated impedance matrix frequency (" << (m->_frequency*1e-6) << " MHz) is more than 2 MHz away from desired frequency (" << (frequency*1e-6) << " MHz).\n";
		std::cout << "Selected impedance table for frequency " << (m->_frequency*1e-6) << " MHz\n";

		return m->_values;
	}
	
	static void Get(double frequency, std::complex<double>* dest)
	{
		memcpy(dest, Get(frequency), sizeof(std::complex<double>)*32*32);
	}
	
private:
	typedef std::complex<double> ctype;
	
	struct ImpedanceMatrix {
		ImpedanceMatrix(double frequency, const ctype* values) :
			_frequency(frequency)
		{
			memcpy(_values, values, sizeof(std::complex<double>)*32*32);
		}
		
		double _frequency;
		ctype _values[32*32];
	};
	
	static void initializeMatrices();
	
	static std::vector<ImpedanceMatrix> matrices;
};

#endif
