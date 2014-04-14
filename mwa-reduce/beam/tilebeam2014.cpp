#include "tilebeam2014.h"

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

TileBeam2014::TileBeam2014(const double* delays) :
	_dipoleHeight(0.28), /* Seems to be 0.3 in the RTS, 0.278 in beam script */
	_dipoleSeparations(1.100),
	_delayStep(435.0e-12)
{
	const double dipoleEast[16] = {-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5,-1.5,-0.5,0.5,1.5};
	const double dipoleNorth[16] = {1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-1.5,-1.5};
	for(size_t i=0;i!=16;++i)
	{
		_dipoleNorth[i] = dipoleNorth[i] * _dipoleSeparations;
		_dipoleEast[i] = dipoleEast[i] * _dipoleSeparations;
		_delays[i] = delays[i]*SPEED_OF_LIGHT*_delayStep;
	}
}

void TileBeam2014::invert32x32(const std::complex<double>* input, std::complex<double>* output)
{
	gsl_matrix_complex
		*gslInput = gsl_matrix_complex_alloc(32, 32),
		*gslOutput = gsl_matrix_complex_alloc(32, 32);
	gsl_permutation *perm = gsl_permutation_alloc(32);
	
	memcpy(gsl_matrix_complex_ptr(gslInput, 0, 0), input, sizeof(std::complex<double>)*32*32);
	
		// Make LU decomposition of matrix m
	int s;
	gsl_linalg_complex_LU_decomp(gslInput, perm, &s);

	// Invert the matrix m
	gsl_linalg_complex_LU_invert(gslInput, perm, gslOutput);
	
	memcpy(output, gsl_matrix_complex_ptr(gslOutput, 0, 0), sizeof(std::complex<double>)*32*32);
	
	gsl_matrix_complex_free(gslInput);
	gsl_matrix_complex_free(gslOutput);
	gsl_permutation_free(perm);
}
