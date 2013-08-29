#include "solutionfile.h"
#include "solutionapplier.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "stdint.h"

using namespace casa;

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std::cout << "Usage: applysolutions [s xx xy yx yy] <ms> <gains-bin-file>\n"
				"Will apply the found solution matrices.\n";
    } else {
		size_t argi = 1;
		double xx=0.0, xy=0.0, yx=0.0, yy=0.0;
		bool preset = false;
		if(strcmp(argv[argi], "-s") == 0)
		{
			xx = atof(argv[argi+1]);
			xy = atof(argv[argi+2]);
			yx = atof(argv[argi+3]);
			yy = atof(argv[argi+4]);
			argi += 5;
			preset = true;
		}
		MeasurementSet ms(argv[argi], Table::Update);
		SolutionFile solutionFile;
		solutionFile.OpenForReading(argv[2]);
		
		SolutionApplier applier;
		if(preset)
			applier.SetPresets(xx, xy, yx, yy);
		applier.Apply(ms, solutionFile);
	}
}
