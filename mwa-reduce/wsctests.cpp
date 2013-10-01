#include "cleanalgorithm.h"
#include "uvector.h"

#include "aocommon/stopwatch.h"

#include <iostream>
#include <random>

int main(int argc, char* argv[])
{
	size_t n = 5120;
	ao::uvector<double> img(n * n), psf(n * n, 0.0);
	std::mt19937 mt;
	mt.seed(42);
	std::normal_distribution<double> normal_dist(0.0, 1.0);
	for(size_t i=0; i!=n*n; ++i)
	{
		img[i] = normal_dist(mt);
	}
	
	size_t x, y;
	volatile size_t copyX, copyY;
	Stopwatch watch(true);
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		CleanAlgorithm::PartialFindPeak(img.data(), n, n, x, y, true, 0, n/2);
		copyX = x; copyY = y;
	}
	std::cout << "PartialFindPeak: " << watch.ToMilliSecondsString() << '\n';
	x=n/2; y=n/2;
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
		CleanAlgorithm::PartialSubtractImage(img.data(), n, n, psf.data(), n, n, x, y, 0.5, 0, n/2);
	std::cout << "PartialSubtractImage: " << watch.ToMilliSecondsString() << '\n';
}
