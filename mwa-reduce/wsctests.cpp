#include "cleanalgorithm.h"
#include "uvector.h"

#include "aocommon/stopwatch.h"

#include <iostream>
#include <random>

void test(bool test, const char* descr)
{
	if(test)
		std::cout << '.' << std::flush;
	else
		std::cout << "\nFAILURE: " << descr << "\n";
}

int main(int argc, char* argv[])
{
	ao::uvector<double> img(16, 0);
	size_t x=size_t(-1), y=size_t(-1);
#ifdef __AVX__
	img[0] = 1;
	CleanAlgorithm::PartialFindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2);
	test(x == 0 && y == 0, "x,y == 0,0");
	img[1] = 2;
	CleanAlgorithm::PartialFindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2);
	test(x == 1 && y == 0, "x,y == 1,0");
	img[4] = 3;
	CleanAlgorithm::PartialFindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2);
	test(x == 0 && y == 1, "x,y == 0,1");
	img[7] = 4;
	CleanAlgorithm::PartialFindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2);
	test(x == 3 && y == 1, "x,y == 3,1");
	img[15] = 5;
	CleanAlgorithm::PartialFindPeakAVX(img.data(), 4, 2, x, y, true, 0, 4);
	test(x == 3 && y == 3, "x,y == 3,3");
	std::cout << '\n';
#endif

	size_t n = 5120;
	ao::uvector<double> psf(n * n, 0.0);
	img.resize(n * n);
	std::mt19937 mt;
	mt.seed(42);
	std::normal_distribution<double> normal_dist(0.0, 1.0);
	for(size_t i=0; i!=n*n; ++i)
	{
		img[i] = normal_dist(mt);
	}
	
	volatile size_t copyX, copyY;
	
	Stopwatch watch(true);
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		CleanAlgorithm::PartialFindPeak(img.data(), n, n, x, y, true, 0, n/2);
		copyX = x; copyY = y;
	}
	std::cout << "PartialFindPeak: " << watch.ToMilliSecondsString() << '\n';
	
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		CleanAlgorithm::PartialFindPeakSimple(img.data(), n, n, x, y, true, 0, n/2);
		copyX = x; copyY = y;
	}
	std::cout << "PartialFindPeakSimple: " << watch.ToMilliSecondsString() << '\n';
	
#ifdef __AVX__
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		CleanAlgorithm::PartialFindPeakAVX(img.data(), n, n, x, y, true, 0, n/2);
		copyX = x; copyY = y;
	}
	std::cout << "PartialFindPeakAVX: " << watch.ToMilliSecondsString() << '\n';
#endif
	
	x=n/2; y=n/2;
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
		CleanAlgorithm::PartialSubtractImage(img.data(), n, n, psf.data(), n, n, x, y, 0.5, 0, n/2);
	std::cout << "PartialSubtractImage: " << watch.ToMilliSecondsString() << '\n';
	
#ifdef __AVX__
	x=n/2; y=n/2;
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
		CleanAlgorithm::PartialSubtractImageAVX(img.data(), n, n, psf.data(), n, n, x, y, 0.5, 0, n/2);
	std::cout << "PartialSubtractImageAVX: " << watch.ToMilliSecondsString() << '\n';
#endif

}
