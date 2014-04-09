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
      std::cout << "Usage: applysolutions [-copy] [-s xx xy yx yy] <ms> <gains-bin-file>\n"
				"Will apply the found solution matrices.\n";
    } else {
		size_t argi = 1;
		double xx=0.0, xy=0.0, yx=0.0, yy=0.0;
		bool preset = false, copyData = false;
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "s")
			{
				xx = atof(argv[argi+1]);
				xy = atof(argv[argi+2]);
				yx = atof(argv[argi+3]);
				yy = atof(argv[argi+4]);
				argi += 4;
				preset = true;
			}
			else if(p == "copy")
			{
				copyData = true;
			}
			++argi;
		}
		MeasurementSet ms(argv[argi], Table::Update);
		SolutionFile solutionFile;
		solutionFile.OpenForReading(argv[argi+1]);
		
		SolutionApplier applier;
		if(preset)
			applier.SetPresets(xx, xy, yx, yy);
		if(copyData)
			applier.SetOutputColumn(casa::MeasurementSet::columnName(casa::MSMainEnums::CORRECTED_DATA));
		applier.Apply(ms, solutionFile);
	}
}
