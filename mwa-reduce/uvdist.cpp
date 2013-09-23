#include <iostream>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>
#include "uvwdistribution.h"

int main(int argc, char* argv[])
{
	if(argc <= 1)
	{
		std::cout << "uvdist: Measure distribution of uvw vector magnitude.\nSyntax:\n\tuvdist <measurementset>";
	}
	else {
		casa::MeasurementSet ms(argv[1]);
		UvwDistribution uvwDist;
		uvwDist.Calculate(ms);
		
		std::cout << "Min uvw dist: " << uvwDist.MinUvw() << ", max uvw dist: " << uvwDist.MaxUvw() << "\n";
		
		std::ofstream file("uvw-dist-data.txt");
		const double minUvw = uvwDist.MinUvw(), maxUvw = uvwDist.MaxUvw();
		const double uvwRange = maxUvw - minUvw;
		for(size_t i=0; i!=uvwDist.BinCount()*10; ++i)
		{
			double distance = (double(i)*0.1 * uvwRange / double(uvwDist.BinCount()) + minUvw);
			//file << distance << '\t' << uvwDist.CountWithInterpolation(distance) << '\n';
			file << distance << '\t' << uvwDist.WeightFromFit(distance) << '\n';
		}
		
		double e, f;
		uvwDist.FitPowerlaw(e, f);
		
		std::ofstream plotFile("uvw-dist.plt");
		plotFile <<
			"set terminal postscript enhanced color\n"
			"set logscale y\n"
			"set xrange [" << minUvw << ':' << maxUvw << "]\n"
			"#set yrange [..:..]\n"
			"set output \"uvw-dist.ps\"\n"
			"set key bottom left\n"
			"set xlabel \"UVW distance (lambda)\"\n"
			"set ylabel \"Count\"\n"
			"plot \"uvw-dist-data.txt\" using 1:2 with lines lw 2.0 title \"\",\\\n"
			"2.718281**((" << f << ")*(x**" << e << ")) with lines lw 2.0 title \"\"\n";
	}
}
