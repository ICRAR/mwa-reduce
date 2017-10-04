#include <iostream>
#include <cmath>
#include <fstream>

#include "mrccatalogue.h"

#include "units/imagecoordinates.h"

#include "model/model.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: mrc2model <mrc-catalogue> <model>\n";
	} else {
		const char *catFilename = argv[1];
		const char *modelFilename = argv[2];
		
		MRCCatalogue catalogue(catFilename);

		Model model;
		ModelSource source;
		while(catalogue.ReadNext(source))
		{
			model.AddSource(source);
		}
		std::cout << "Saving model with " << model.SourceCount() << " sources.\n";
		model.Save(modelFilename);
	}
}
