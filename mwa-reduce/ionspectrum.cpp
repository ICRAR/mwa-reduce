#include "ionspectrummaker.h"

int main(int argc, char* argv[])
{
	if(argc < 4) {
		std::cout << "syntax: ionspectrum <ms> <ionsolutions> <model>\n";
		return -1;
	}
	
	const char *msFilename(argv[1]);
	const char *ionFilename(argv[2]);
	const char *modelFilename(argv[3]);
	IonSpectrumMaker isMaker;
	isMaker.Accumulate(msFilename, ionFilename, modelFilename);
	
	return 0;
}
