#include "ionspectrummaker.h"

int main(int argc, char* argv[])
{
	if(argc < 4) {
		std::cout << "syntax: ionspectrum <ms> <ionsolutions> <model> <output>\n";
		return -1;
	}
	
	const char *msFilename(argv[1]);
	const char *ionFilename(argv[2]);
	const char *modelFilename(argv[3]);
	const char *outputFilename(argv[4]);
	IonSpectrumMaker isMaker(msFilename, ionFilename, modelFilename);
	isMaker.Accumulate();
	isMaker.Save(outputFilename);
	
	return 0;
}
