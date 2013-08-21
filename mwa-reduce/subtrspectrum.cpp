#include <iostream>

#include "spectrumsubtractor.h"

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout << "Syntax: subtrspectrum <model> <ms>\nDetermines spectrum for each source in the model at each time step,\nand subtracts it from the measurement set.";
		return -1;
	}
	
	Model model(argv[1]);
	casa::MeasurementSet ms(argv[2], casa::MeasurementSet::Update);
	SpectrumSubtractor subtractor(ms, model);
	subtractor.Perform();
}
