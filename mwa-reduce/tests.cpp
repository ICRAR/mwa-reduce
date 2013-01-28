#include <iostream>

#include "sourcesdfwithsamples.h"

int main(int argc, char *argv[])
{
	SourceSDFWithSamples<long double> sdf1;
	sdf1.AddSample(1.0, 2.0);
	sdf1.AddSample(1.0, 3.0);
	std::cout << "SDF1@1 : expecting 1, is " << sdf1.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF1@2 : expecting 1, is " << sdf1.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF1@3 : expecting 1, is " << sdf1.FluxAtFrequency(3.0) << '\n';
	std::cout << "SDF1@4 : expecting 1, is " << sdf1.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF1@2.5 : expecting 1, is " << sdf1.FluxAtFrequency(2.5) << '\n';
	
	std::cout << "SDF1 int[2,3] : expecting 1, is " << sdf1.IntegratedFlux(2.0, 3.0) << '\n';
	std::cout << "SDF1 int[1,3] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 3.0) << '\n';
	std::cout << "SDF1 int[1,2] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF1 int[3,4] : expecting 1, is " << sdf1.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF1 int[2,4] : expecting 1, is " << sdf1.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF1 int[1,4] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF1 int[1,1] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 1.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf2;
	sdf2.AddSample(1.0, 2.0);
	sdf2.AddSample(2.0, 4.0);
	std::cout << "SDF2@1 : expecting 0.5, is " << sdf2.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF2@2 : expecting 1, is " << sdf2.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF2@4 : expecting 2, is " << sdf2.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF2@8 : expecting 4, is " << sdf2.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF2@3 : expecting 1.5, is " << sdf2.FluxAtFrequency(3.0) << '\n';
	
	std::cout << "SDF2 int[2,4] : expecting 1.5, is " << sdf2.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,2] : expecting 0.75, is " << sdf2.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF2 int[3,4] : expecting 1.75, is " << sdf2.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,4] : expecting 1.25, is " << sdf2.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,1] : expecting 0.5, is " << sdf2.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF2 int[3,5] : expecting 2, is " << sdf2.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF2 int[1,8] : expecting 2.25, is " << sdf2.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf3(sdf2);
	sdf3.AddSample(2.0, 6.0);
	std::cout << "SDF3@1 : expecting 0.5, is " << sdf3.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF3@2 : expecting 1, is " << sdf3.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF3@4 : expecting 2, is " << sdf3.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF3@6 : expecting 2, is " << sdf3.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF3@8 : expecting 2, is " << sdf3.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF3@3 : expecting 1.5, is " << sdf3.FluxAtFrequency(3.0) << '\n';

	std::cout << "SDF3 int[2,4] : expecting 1.5, is " << sdf3.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,2] : expecting 0.75, is " << sdf3.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF3 int[3,4] : expecting 1.75, is " << sdf3.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,4] : expecting 1.25, is " << sdf3.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,1] : expecting 0.5, is " << sdf3.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF3 int[3,5] : expecting 1.875, is " << sdf3.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF3 int[1,8] : expecting 1.679, is " << sdf3.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf4;
	sdf4.AddSample(-1.0, 2.0);
	sdf4.AddSample(-2.0, 4.0);
	sdf4.AddSample(-2.0, 6.0);
	std::cout << "SDF4@1 : expecting -0.5, is " << sdf4.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF4@2 : expecting -1, is " << sdf4.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF4@4 : expecting -2, is " << sdf4.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF4@6 : expecting -2, is " << sdf4.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF4@8 : expecting -2, is " << sdf4.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF4@3 : expecting -1.5, is " << sdf4.FluxAtFrequency(3.0) << '\n';
	
	std::cout << "SDF4 int[2,4] : expecting -1.5, is " << sdf4.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,2] : expecting -0.75, is " << sdf4.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF4 int[3,4] : expecting -1.75, is " << sdf4.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,4] : expecting -1.25, is " << sdf4.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,1] : expecting -0.5, is " << sdf4.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF4 int[3,5] : expecting -1.875, is " << sdf4.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF4 int[1,8] : expecting -1.679, is " << sdf4.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf5;
	sdf5.AddSample(2.0, 2.0);
	sdf5.AddSample(-2.0, 4.0);
	sdf5.AddSample(-2.0, 6.0);
	sdf5.AddSample(2.0, 8.0);
	std::cout << "SDF5@1 : expecting 4, is " << sdf5.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF5@2 : expecting 2, is " << sdf5.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF5@3 : expecting 0, is " << sdf5.FluxAtFrequency(3.0) << '\n';
	std::cout << "SDF5@4 : expecting -2, is " << sdf5.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF5@6 : expecting -2, is " << sdf5.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF5@7 : expecting 0, is " << sdf5.FluxAtFrequency(7.0) << '\n';
	std::cout << "SDF5@8 : expecting 2, is " << sdf5.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF5@9 : expecting 4, is " << sdf5.FluxAtFrequency(9.0) << '\n';
	
	std::cout << "SDF5 int[2,4] : expecting 0, is " << sdf5.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF5 int[0,2] : expecting 4, is " << sdf5.IntegratedFlux(0.0, 2.0) << '\n';
	std::cout << "SDF5 int[0,4] : expecting 2, is " << sdf5.IntegratedFlux(0.0, 4.0) << '\n';
}
