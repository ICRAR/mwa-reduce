#include <iostream>
#include <stdexcept>
#include "peeler.h"

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: peel [-datacolumn <column>] [-beam-on-source] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist>] [-l <precision>] [-i <niter>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-applybeam] [-t <solution interval timesteps>] <measurementset.ms>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n";
	} else {
		int argi = 1;
		bool
			beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false;
		std::string modelFile, rhsSolutionFile;
		std::string dataColumnName = "DATA";
		size_t niter = 25, solutionInterval = 1;
		double limit = 0.0001, minUVW = 0.0;
		
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-i") == 0)
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-l") == 0)
			{
				limit = atof(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-m") == 0)
			{
				modelFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-minuv") == 0)
			{
				minUVW = atof(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-applybeam") == 0)
			{
				applyBeam = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-beam-on-source") == 0)
			{
				beamOnSource = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-scalar") == 0)
			{
				onlyScalar = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-diag") == 0)
			{
				onlyDiag = true;
				++argi;
			}
			else if(strcmp(argv[argi], "-rotation") == 0)
			{
				onlyRotation = true;
				argi++;
			}
			else if(strcmp(argv[argi], "-rhs") == 0)
			{
				rhsSolutionFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-datacolumn") == 0)
			{
				dataColumnName = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-t") == 0)
			{
				solutionInterval = atoi(argv[argi+1]);
				argi += 2;
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
		}
		
		if(argc <= argi) throw std::runtime_error("Incorrect number of parameters");
		
		const char *msName = argv[argi];
		//const char *outName = argv[argi+1];
		
		casa::MeasurementSet ms(msName);
		Peeler peeler(ms);
		
		peeler.SetNIter(niter);
		peeler.SetLimit(limit);
		peeler.SetModelFilename(modelFile);
		peeler.SetMinUVW(minUVW);
		peeler.SetApplyBeam(applyBeam);
		peeler.SetBeamOnSource(beamOnSource);
		peeler.SetOnlyScalar(onlyScalar);
		peeler.SetOnlyDiag(onlyDiag);
		peeler.SetOnlyRotation(onlyRotation);
		peeler.SetRHSSolutionFile(rhsSolutionFile);
		peeler.SetDataColumnName(dataColumnName);
		peeler.SetSolutionInterval(solutionInterval);
		
		peeler.Perform();
	}
}
