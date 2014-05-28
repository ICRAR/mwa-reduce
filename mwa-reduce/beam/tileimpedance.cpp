#include "tileimpedance.h"

#include <fitsio.h>

std::vector<TileImpedance::ImpedanceMatrix> TileImpedance::matrices = std::vector<TileImpedance::ImpedanceMatrix>();

void TileImpedance::initializeMatrices()
{
	int status = 0;
	fitsfile *fitsPtr;
	fits_open_file(&fitsPtr, "ZMatrix.fits", READONLY, &status);
	
	int hdunum = 0;
	fits_get_num_hdus(fitsPtr, &hdunum, &status);
	
	std::cout.precision(16);
	std::complex<double> cValues[32*32];
	for(int hdu=1; hdu<=hdunum; ++hdu)
	{
		fits_movabs_hdu(fitsPtr, hdu, NULL, &status);
		
		char keyValue[FLEN_VALUE];
		fits_read_keyword(fitsPtr, "FREQ", keyValue, NULL, &status);
		double frequency = atof(keyValue);
		
		double values[2*32*32];
		
		long fpixel[3] = {1, 1, 1};
		fits_read_pix(fitsPtr, TDOUBLE, fpixel, 2*32*32, 0, values, 0, &status);
		
		for(size_t v=0; v!=32*32; ++v)
		{
			double mag=values[v], ph = values[v + 32*32];
			double real = mag*cos(ph), imag = mag*sin(ph);
			cValues[v] = std::complex<double>(real, imag);
		}
		ImpedanceMatrix matrix(frequency, cValues);
		matrices.push_back(matrix);
	}
}
