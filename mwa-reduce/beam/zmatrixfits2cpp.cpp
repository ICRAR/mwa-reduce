#include <fitsio.h>
#include <iostream>
#include <limits>
#include <cmath>

int main(int argc, char* argv[])
{
	int status = 0;
	fitsfile *fitsPtr;
	fits_open_file(&fitsPtr, argv[1], READONLY, &status);
	
	int hdunum = 0;
	fits_get_num_hdus(fitsPtr, &hdunum, &status);
	
	std::cout.precision(16);
	for(int hdu=1; hdu<=hdunum; ++hdu)
	{
		fits_movabs_hdu(fitsPtr, hdu, NULL, &status);
		
		char keyValue[FLEN_VALUE];
		fits_read_keyword(fitsPtr, "FREQ", keyValue, NULL, &status);
		double frequency = atof(keyValue);
		
		double values[2*32*32];
		
		int anynul = 0;
		long fpixel[3] = {1, 1, 1};
		fits_read_pix(fitsPtr, TDOUBLE, fpixel, 2*32*32, 0, values, 0, &status);
		
		std::cout << "ImpedanceMatrix(" << frequency << ",\n{\n  ";
		for(size_t v=0; v!=32*32; ++v)
		{
			double mag=values[v], ph = values[v + 32*32];
			double real = mag*cos(ph), imag = mag*sin(ph);
			if(v != 0)
			{
				std::cout << ',';
				if(v%2 == 0) std::cout << "\n  ";
				else std::cout << ' ';
			}
			std::cout << "ctype(" << real << ", " << imag << ")";
		}
		std::cout << "\n})";
		if(hdu != hdunum) std::cout << ',';
		std::cout << "\n";
	}
}
