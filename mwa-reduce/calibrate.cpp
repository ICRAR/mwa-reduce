#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>
#include "calibrator.h"
#include "calibrationmethod.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout
			<< "Usage: calibrate [-beam-on-source] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist in m>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-j <threads>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-applybeam] [-t timesteps] [-datacolumn <name>] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n";
	} else {
		int argi = 1;
		// bool saveCrossTermsPlotFile = false, saveFaradayPlotFiles = false;
		bool
			savePlotFiles = false, beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, modelFile, rhsSolutionFile;
		size_t niter = CalibrationMethod::DefaultNIter(), solutionInterval = 0;
		std::string dataColumnName = "DATA";
		double
			minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
			stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy(),
			minUVW = 0.0;
		size_t threadCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		
		while(argv[argi][0] == '-')
		{
			std::string param(&argv[argi][1]);
			if(param == "p")
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 3;
			}
			else if(param == "datacolumn")
			{
				dataColumnName = argv[argi+1];
				argi += 2;
			}
			/*else if(param == "pf") == 0)
			{
				saveFaradayPlotFiles = true;
				plotFaradayFile = argv[argi+1];
				argi += 2;
			}
			else if(param == "px") == 0)
			{
				saveCrossTermsPlotFile = true;
				crossTermsPlotFile = argv[argi+1];
				argi += 2;
			}*/
			else if(param == "i")
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "j")
			{
				threadCount = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "a")
			{
				minAccuracy = atof(argv[argi+1]);
				stopAccuracy = atof(argv[argi+2]);
				argi += 3;
			}
			else if(param == "m")
			{
				modelFile = argv[argi+1];
				argi += 2;
			}
			else if(param == "t")
			{
				solutionInterval = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(param == "minuv")
			{
				minUVW = atof(argv[argi+1]);
				argi += 2;
			}
			else if(param == "applybeam")
			{
				applyBeam = true;
				++argi;
			}
			else if(param == "beam-on-source")
			{
				beamOnSource = true;
				++argi;
			}
			else if(param == "scalar")
			{
				onlyScalar = true;
				++argi;
			}
			else if(param == "diag")
			{
				onlyDiag = true;
				++argi;
			}
			else if(param == "rhs")
			{
				rhsSolutionFile = argv[argi+1];
				argi += 2;
			}
			else if(param == "rotation")
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
		
		Calibrator calibrator(ms, threadCount);
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
		calibrator.SetSavePlotFiles(savePlotFiles);
		calibrator.SetDataColumnName(dataColumnName);
		if(savePlotFiles)
		{
			calibrator.SetPlotFilenames(plotPhaseFile, plotGainFile);
		}
		calibrator.Perform();
	}
}
