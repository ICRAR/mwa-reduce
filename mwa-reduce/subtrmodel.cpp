#include <iostream>
#include <stdexcept>

#include <ms/MeasurementSets/MeasurementSet.h>

#include "model.h"
#include "subtractor.h"

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: subtrmodel [-usemodelcol] [-datacolumn <COLUMN>] [-applybeam] [-r / -s] [-n <σ>] <model> <ms>\n"
			"Subtracts the model from the visibilities. This 'peels' the\n"
			"sources out. Only affects cross-correlations. -r to revert or -s to set.\n";
	} else {
		bool revert = false , setvis = false, addNoise = false, applyBeam = false, useModelCol = false;
		double noiseSigma = 1.0;
		size_t argi = 1;
		size_t threadCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		std::string dataColumnName = "DATA";
		
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-r") == 0) { revert=true; }
			else if(strcmp(argv[argi], "-s") == 0) { setvis=true; }
			else if(strcmp(argv[argi], "-n") == 0) { addNoise=true; ++argi; noiseSigma = atof(argv[argi]); }
			else if(strcmp(argv[argi], "-applybeam") == 0) { applyBeam=true; }
			else if(strcmp(argv[argi], "-datacolumn") == 0) { ++argi; dataColumnName=argv[argi]; }
			else if(strcmp(argv[argi], "-usemodelcol") == 0) { useModelCol=true; }
			else throw std::runtime_error("Invalid param");
			++argi;
		}
		
		std::cout << "Reading model... " << std::flush;
		Model model(argv[argi]);
		std::cout << "DONE\n";
		
		std::cout << "Opening measurement set... " << std::flush;
		MeasurementSet ms(argv[argi+1], Table::Update);
		std::cout << "DONE\n";
		
		if(useModelCol)
		{
			ProgressBar progress("Subtracting model column");
			casa::ArrayColumn<casa::Complex> dataColumn(ms, dataColumnName);
			casa::ROArrayColumn<casa::Complex> modelColumn(ms, casa::MeasurementSet::columnName(casa::MeasurementSet::MODEL_DATA));
			casa::Array<casa::Complex>
				dataArr(dataColumn.shape(0)), modelArr(modelColumn.shape(0));
			for(size_t i=0; i!=ms.nrow(); ++i)
			{
				dataColumn.get(i, dataArr);
				modelColumn.get(i, modelArr);
				
				casa::Array<casa::Complex>::contiter d=dataArr.cbegin();
				for(casa::Array<casa::Complex>::const_contiter m=modelArr.cbegin(); m!=modelArr.cend(); ++m)
				{
					*d -= *m;
				}
				dataColumn.put(i, dataArr);
				progress.SetProgress(i+1, ms.nrow());
			}
		}
		else {
			Subtractor subtractor(threadCount);
			subtractor.SetRevert(revert);
			subtractor.SetToModel(setvis);
			subtractor.SetAddNoise(addNoise);
			subtractor.SetApplyBeam(applyBeam);
			subtractor.SetNoiseSigma(noiseSigma);
			subtractor.SetDataColumn(dataColumnName);
			subtractor.Subtract(ms, model);
		}
	}
}
