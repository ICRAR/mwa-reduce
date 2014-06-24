
#include "tileimpedance.h"
#include "lnaimpedance.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

static_assert(sizeof(std::complex<double>) == sizeof(gsl_complex), "sizeof(std::complex<double>) != sizeof(gsl_complex)");

int main(int argc, char* argv[])
{
	double freqs[]={80,130,150,200,230};
	double delays[] = {6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0,6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0};
	double dq=435e-12*3e8;
	double za=14;
	for(size_t i=0; i!=5; ++i)
	{
		double fs = freqs[i]*1e6;
		std::complex<double> lnaImp = LNAImpedance::Get(fs);
    std::cout << "LNA Impedance at " << freqs[i] << " MHz: " <<  lnaImp << '\n';
    
		double lam = 3e8/fs;
		gsl_matrix_complex *impMatrix = gsl_matrix_complex_alloc(32, 32);
		std::complex<double> *impMatrixPtr = reinterpret_cast<std::complex<double>*>(gsl_matrix_complex_ptr(impMatrix, 0, 0));
		
		TileImpedance::Get(fs, impMatrixPtr);
		std::complex<double> ph_rot[32];
		for(size_t t=0; t!=32; ++t)
		{
			double phase = M_PI*-2.0*delays[t]*(dq/lam);
			ph_rot[t] = std::complex<double>(cos(phase), sin(phase));
		}
		
		// Add lna impedance to diagonal values
		for(size_t diag=0; diag!=32; ++diag)
			impMatrixPtr[diag*33] += lnaImp;
		
		std::cout << "MWA " << freqs[i] << + " MHz total Z matrix magnitude\n";
		for(size_t y=0; y!=32; ++y) {
			for(size_t x=0; x!=32; ++x) {
				std::complex<double> val = *reinterpret_cast<std::complex<double>*>(gsl_matrix_complex_ptr(impMatrix, y, x));
				std::cout << round(std::abs(val)*10.0)/10.0 << ' ';
			}
			std::cout << '\n';
		}
		
		gsl_matrix_complex *inverse = gsl_matrix_complex_alloc(32, 32);
		gsl_permutation *perm = gsl_permutation_alloc(32);
		
			// Make LU decomposition of matrix m
		int s;
		gsl_linalg_complex_LU_decomp(impMatrix, perm, &s);

		// Invert the matrix m
		gsl_linalg_complex_LU_invert(impMatrix, perm, inverse);
		
		
		std::complex<double> current[32];
		for(size_t j=0; j!=32; ++j)
		{
			for(size_t i=0; i!=32; ++i)
				current[j] += *reinterpret_cast<std::complex<double>*>(gsl_matrix_complex_ptr(inverse, j, i)) * ph_rot[i];
		}
		std::cout << "MWA " << freqs[i] << "MHz X dipole current amplitude (ZA=" << za << "deg)\n";
		for(size_t y=0; y!=4; ++y) {
			for(size_t x=0; x!=4; ++x) {
				std::cout << std::abs(current[16+y*4 + x])*1000.0 << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "MWA " << freqs[i] << "MHz Y dipole current amplitude (ZA=" << za << "deg)\n";
		for(size_t y=0; y!=4; ++y) {
			for(size_t x=0; x!=4; ++x) {
				std::cout << std::abs(current[y*4 + x])*1000.0 << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "MWA " << freqs[i] << "MHz X dipole current phase (ZA=" << za << "deg)\n";
		for(size_t y=0; y!=4; ++y) {
			for(size_t x=0; x!=4; ++x) {
				std::cout << std::arg(current[16+y*4 + x])*180.0/M_PI << ' ';
			}
			std::cout << '\n';
		}
		std::cout << "MWA " << freqs[i] << "MHz Y dipole current phase (ZA=" << za << "deg)\n";
		for(size_t y=0; y!=4; ++y) {
			for(size_t x=0; x!=4; ++x) {
				std::cout << std::arg(current[y*4 + x])*180.0/M_PI << ' ';
			}
			std::cout << '\n';
		}
		
		//current = numpy.dot(inv_z,ph_rot).reshape(2,4,4)
		gsl_matrix_complex_free(inverse);
		gsl_matrix_complex_free(impMatrix);
		gsl_permutation_free(perm);
	}
}
