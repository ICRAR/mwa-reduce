#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>
#include "calibrator.h"
#include "calibrationmethod.h"

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: calibrate [-beam-on-source] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-applybeam] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n";
	} else {
		int argi = 1;
		bool
			savePlotFiles = false, saveFaradayPlotFiles = false, saveCrossTermsPlotFile = false, beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, modelFile, rhsSolutionFile;
		size_t niter = CalibrationMethod::DefaultNIter(), solutionInterval = 0;
		double
			minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
			stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy(),
			minUVW = 0.0;
		
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-p") == 0)
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 3;
			}
			else if(strcmp(argv[argi], "-pf") == 0)
			{
				saveFaradayPlotFiles = true;
				plotFaradayFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-px") == 0)
			{
				saveCrossTermsPlotFile = true;
				crossTermsPlotFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-i") == 0)
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-a") == 0)
			{
				minAccuracy = atof(argv[argi+1]);
				stopAccuracy = atof(argv[argi+2]);
				argi += 3;
			}
			else if(strcmp(argv[argi], "-m") == 0)
			{
				modelFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-t") == 0)
			{
				solutionInterval = atoi(argv[argi+1]);
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
			else if(strcmp(argv[argi], "-rhs") == 0)
			{
				rhsSolutionFile = argv[argi+1];
				argi += 2;
			}
			else if(strcmp(argv[argi], "-rotation") == 0)
			{
				onlyRotation = true;
				argi++;
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
		}
		
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");
		
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		casa::MeasurementSet ms(msName);
		
		Calibrator calibrator(ms);
		calibrator.SetNIter(niter);
		calibrator.SetAccuracy(minAccuracy, stopAccuracy);
		calibrator.SetModelFilename(modelFile);
		calibrator.SetSolutionInterval(solutionInterval);
		calibrator.SetMinUVW(minUVW);
		calibrator.SetApplyBeam(applyBeam);
		calibrator.SetBeamOnSource(beamOnSource);
		calibrator.SetOnlyScalar(onlyScalar);
		calibrator.SetOnlyDiag(onlyDiag);
		calibrator.SetRHSSolutionFile(rhsSolutionFile);
		calibrator.SetOnlyRotation(onlyRotation);
		calibrator.SetSolutionOutputFilename(outName);
		calibrator.SetVerbose(true);
		calibrator.Perform();
	}
}
