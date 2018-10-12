#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>

#include "calibrator.h"
#include "calibrationmethod.h"
#include "system.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Usage: calibrate [-absmem <mem in GB>] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist in m>] [-maxuv <max uvw dist in m>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-j <threads>] [-m <model>] [-scalar] [-diag] [-rotation] [-applybeam] [-t timesteps] [-interval <starttime> <endtime>] [-ch channels per solution] [-datacolumn <name>] [-quiet] [-mwa-path path] <measurementset.ms> <solutions.bin>\n\n"
			"This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n\n"
			"The official name of this algorithm is the \"Mitchcal\" algorithm. The following is a suggestion for referencing this algorithm in scientific articles:\n"
			"\" Calibration was performed with the full-Jones Mitchcal algorithm developed for MWA calibration, as described by Offringa et al. (2016) \"\n"
			"bibtex:\n"
			"@article{offringa-2016,\n"
			"  author = {Offringa, A. R. and Trott, C. M. and Hurley-Walker and others},\n"
			"  title = {Parametrizing Epoch of Reionization foregrounds: a deep survey of low-frequency point-source spectra with the Murchison Widefield Array},\n"
			"  volume = {458}, \n"
			"  number = {1}, \n"
			"  pages = {1057-1070}, \n"
			"  year = {2016}, \n"
			"  doi = {10.1093/mnras/stw310}, \n"
			"  URL = {http://mnras.oxfordjournals.org/content/458/1/1057.abstract}, \n"
			"  journal = {MNRAS}\n}\n";
	} else {
		int argi = 1;
		// bool saveCrossTermsPlotFile = false, saveFaradayPlotFiles = false;
		bool
			savePlotFiles = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false, doQuiet = false,
			hasInterval = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, modelFile,  mwaPath;
		size_t
			niter = CalibrationMethod::DefaultNIter(),
			solutionInterval = 0,
			solutionChannels = 1;
		std::string dataColumnName = "DATA";
		double
			minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
			stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy(),
			minUVW = 0.0,
		  maxUVW = 10000000.0,
		  absmem = 0.0;
		size_t threadCount = System::ProcessorCount();
		size_t intervalStart = 0, intervalEnd = 0;
		size_t followAntenna = 0;
		
		while(argv[argi][0] == '-')
		{
			std::string param(&argv[argi][1]);
			if(param == "p")
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 2;
			}
			else if(param == "datacolumn")
			{
				++argi;
				dataColumnName = argv[argi];
			}
			else if(param == "absmem")
			{
				++argi;
				absmem = atof(argv[argi]);
			}
			else if(param == "i")
			{
				++argi;
				niter = atoi(argv[argi]);
			}
			else if(param == "j")
			{
				++argi;
				threadCount = atoi(argv[argi]);
			}
			else if(param == "a")
			{
				minAccuracy = atof(argv[argi+1]);
				stopAccuracy = atof(argv[argi+2]);
				argi += 2;
			}
			else if(param == "m")
			{
				argi++;
				modelFile = argv[argi];
			}
			else if(param == "t")
			{
				argi++;
				solutionInterval = atoi(argv[argi]);
			}
			else if(param == "ch")
			{
				++argi;
				solutionChannels = atoi(argv[argi]);
			}
			else if(param == "minuv")
			{
				argi++;
				minUVW = atof(argv[argi]);
			}
			else if(param == "maxuv")
			{
				++argi;
				maxUVW = atof(argv[argi]);
			}
			else if(param == "applybeam")
			{
				applyBeam = true;
			}
			else if(param == "scalar")
			{
				onlyScalar = true;
			}
			else if(param == "diag")
			{
				onlyDiag = true;
			}
			else if(param == "rotation")
			{
				onlyRotation = true;
			}
			else if(param == "quiet")
			{
				doQuiet = true;
			}
			else if(param == "follow-antenna")
			{
				++argi;
				followAntenna = atoi(argv[argi]);
			}
			else if(param == "mwa-path")
			{
				++argi;
				mwaPath = argv[argi];
			}
			else if(param == "interval")
			{
				hasInterval = true;
				intervalStart = atof(argv[argi+1]);
				intervalEnd = atof(argv[argi+2]);
				argi += 2;
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
			++argi;
		}
		
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");
		
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		casacore::MeasurementSet ms(msName);
		
		Calibrator calibrator(ms, threadCount);
		calibrator.SetNIter(niter);
		calibrator.SetAccuracy(minAccuracy, stopAccuracy);
		calibrator.SetModelFilename(modelFile);
		calibrator.SetSolutionInterval(solutionInterval);
		calibrator.SetSolutionChannels(solutionChannels);
		calibrator.SetMinUVW(minUVW);
		calibrator.SetMaxUVW(maxUVW);
		calibrator.SetApplyBeam(applyBeam);
		calibrator.SetOnlyScalar(onlyScalar);
		calibrator.SetOnlyDiag(onlyDiag);
		calibrator.SetOnlyRotation(onlyRotation);
		calibrator.SetSolutionOutputFilename(outName);
		calibrator.SetVerbose(!doQuiet);
		calibrator.SetSavePlotFiles(savePlotFiles);
		calibrator.SetDataColumnName(dataColumnName);
		calibrator.SetFollowAntenna(followAntenna);
		calibrator.SetAbsMem(absmem);
		calibrator.SetMWAPath(mwaPath);
		if(hasInterval)
			calibrator.SetInterval(intervalStart, intervalEnd);
		if(savePlotFiles)
		{
			calibrator.SetPlotFilenames(plotPhaseFile, plotGainFile);
		}
		calibrator.Perform();
	}
}
