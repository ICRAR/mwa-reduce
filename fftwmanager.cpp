#include "fftwmanager.h"

#include "system.h"

#include <iostream>

#include <fftw3.h>

FFTWManager::FFTWManager(bool verbose) :
	_multiThreadEnabledDepth(0),
	_verbose(verbose),
	_nThreads(std::min(System::ProcessorCount(), 4u))
{ }

FFTWManager::FFTWManager(size_t nThreads, bool verbose) :
	_multiThreadEnabledDepth(0),
	_verbose(verbose),
	_nThreads(nThreads)
{ }

FFTWManager::~FFTWManager()
{
	fftw_cleanup_threads();
}

void FFTWManager::activateMultipleThreads()
{
	if(_verbose)
		std::cout << "Setting FFTW to use " << _nThreads << " threads.\n";
	fftw_init_threads();
	fftw_plan_with_nthreads(_nThreads);
}

void FFTWManager::endMultipleThreads()
{
	fftw_plan_with_nthreads(1);
}
