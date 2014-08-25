#include "ionspectrummaker.h"

#include <ms/MeasurementSets/MeasurementSet.h>

int main(int argc, char* argv[])
{
	if(argc < 4) {
		std::cout << "syntax: ionspectrumcombine <template-ms> <model> <output prefix> <spectrum1> [<spectrum2> ..]\n";
		return -1;
	}
	casa::MeasurementSet ms(argv[1]);
	const char *modelFilename(argv[2]);

	size_t threadCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
	IonSpectrumMaker isMaker(threadCount);
	BandData bandData(ms.spectralWindow());
	isMaker.InitializeForFileAcc(modelFilename, bandData);
	std::string outputPrefix(argv[3]);
	
	for(int argi=4; argi!=argc; ++argi)
		isMaker.AccumulateFile(argv[argi]);
	
	isMaker.Save((outputPrefix + ".bin").c_str());
	
	Model model;
	isMaker.GetModelWithSpectra(model);
	model.Save(outputPrefix + "-model.txt");
	
	std::cout << "Saved spectra for " << model.SourceCount() << " / " << isMaker.PositionsModel().SourceCount() << " sources.\n";
	
	return 0;
}
