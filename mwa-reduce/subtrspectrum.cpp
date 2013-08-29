#include <iostream>

#include "spectrumsubtractor.h"

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout << "Syntax: subtrspectrum [-datacolumn <column>] [-t <interval>] <model> <ms>\nDetermines spectrum for each source in the model at each time step,\nand subtracts it from the measurement set.";
		return -1;
	}
	int argi = 1;
	std::string dataColumn = "DATA";
	size_t fittingInterval = 1;
	while(argv[argi][0] == '-')
	{
		std::string param(&argv[argi][1]);
		if(param == "datacolumn")
		{
			++argi;
			dataColumn = argv[argi];
		}
		else if(param == "t")
		{
			++argi;
			fittingInterval = atoi(argv[argi]);
		}
		else
			throw std::runtime_error(std::string("Bad param ") + param);
		++argi;
	}
	Model model(argv[argi]);
	casa::MeasurementSet ms(argv[argi+1], casa::MeasurementSet::Update);
	SpectrumSubtractor subtractor(ms, model);
	subtractor.SetDataColumn(dataColumn);
	subtractor.SetFittingInterval(fittingInterval);
	subtractor.Perform();
}
