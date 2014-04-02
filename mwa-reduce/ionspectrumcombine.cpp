#include "ionspectrummaker.h"

int main(int argc, char* argv[])
{
	if(argc < 5) {
		std::cout << "syntax: ionspectrumcombine <ms> <ionsolutions> <model> <output prefix> <spectrum1> [<spectrum2> ..]\n";
		return -1;
	}
	
	const char *msFilename(argv[1]);
	const char *ionFilename(argv[2]);
	const char *modelFilename(argv[3]);
	IonSpectrumMaker isMaker(msFilename, ionFilename, modelFilename);
	std::string outputPrefix(argv[4]);
	
	for(int argi=5; argi!=argc; ++argi)
		isMaker.AccumulateFile(argv[argi]);
	
	isMaker.Save((outputPrefix + ".bin").c_str());
	
	Model model;
	isMaker.GetModel(model);
	model.Save(outputPrefix + "-model.txt");
	
	return 0;
}
