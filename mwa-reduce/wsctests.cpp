#include "cleanalgorithms/simpleclean.h"

#include "uvector.h"

#include "aocommon/stopwatch.h"

#include <iostream>
#include <random>

volatile size_t copyX, copyY;

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
#if defined __AVX__ && !defined FORCE_NON_AVX
	img[0] = 1;
	SimpleClean::FindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2, 0.0);
	test(x == 0 && y == 0, "x,y == 0,0");
	img[1] = 2;
	SimpleClean::FindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2, 0.0);
	test(x == 1 && y == 0, "x,y == 1,0");
	img[4] = 3;
	SimpleClean::FindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2, 0.0);
	test(x == 0 && y == 1, "x,y == 0,1");
	img[7] = 4;
	SimpleClean::FindPeakAVX(img.data(), 4, 2, x, y, true, 0, 2, 0.0);
	test(x == 3 && y == 1, "x,y == 3,1");
	img[15] = 6;
	SimpleClean::FindPeakAVX(img.data(), 4, 4, x, y, true, 0, 4, 0.0);
	test(x == 3 && y == 3, "x,y == 3,3");
	img[14] = 5;
	SimpleClean::FindPeakAVX(img.data(), 3, 5, x, y, true, 0, 5, 0.0);
	test(x == 2 && y == 4, "x,y == 2,4");
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
	
	Stopwatch watch(true);
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		SimpleClean::FindPeak(img.data(), n, n, x, y, true, 0, n/2, 0.0);
		copyX = x; copyY = y;
	}
	std::cout << "FindPeak: " << watch.ToMilliSecondsString() << '\n';
	
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		SimpleClean::FindPeakSimple(img.data(), n, n, x, y, true, 0, n/2, 0.0);
		copyX = x; copyY = y;
	}
	std::cout << "FindPeakSimple: " << watch.ToMilliSecondsString() << '\n';
	
#if defined __AVX__ && !defined FORCE_NON_AVX
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
	{
		SimpleClean::FindPeakAVX(img.data(), n, n, x, y, true, 0, n/2, 0.0);
		copyX = x; copyY = y;
	}
	std::cout << "FindPeakAVX: " << watch.ToMilliSecondsString() << '\n';
#endif
	
	x=n/2; y=n/2;
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
		SimpleClean::PartialSubtractImage(img.data(), n, n, psf.data(), n, n, x, y, 0.5, 0, n/2);
	std::cout << "PartialSubtractImage: " << watch.ToMilliSecondsString() << '\n';
	
#if defined __AVX__ && !defined FORCE_NON_AVX
	x=n/2; y=n/2;
	watch.Reset();
	watch.Start();
	for(size_t repeat=0; repeat!=100; ++repeat)
		SimpleClean::PartialSubtractImageAVX(img.data(), n, n, psf.data(), n, n, x, y, 0.5, 0, n/2);
	std::cout << "PartialSubtractImageAVX: " << watch.ToMilliSecondsString() << '\n';
#endif

}
