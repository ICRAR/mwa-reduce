#include <iostream>
#include "model.h"
#include "spectrummaker.h"

int main(int argc, char* argv[])
{
	if(argc < 4)
		std::cout << "Syntax: combinespectra <positionsmodel> <outputmodel> <intermediate1> [<intermediate2>..]\n";
	else {
		SpectrumMaker spectrumMaker;
		
		Model model(argv[1]);
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
			spectrumMaker.AddSource(*s);
		
		for(size_t argi=3; argi!=argc; ++argi)
		{
			std::cout << "Adding " << argv[argi] << "...\n";
			spectrumMaker.AddMeasurementsFromFile(argv[argi]);
		}
		
		spectrumMaker.Finish();
		
		Model outputModel;
		spectrumMaker.ToModel(outputModel);
		outputModel.Save(argv[2]);
	}
}
