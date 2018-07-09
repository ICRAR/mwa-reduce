#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <iostream>
#include <stdexcept>
#include "calibrator.h"
#include "calibrationmethod.h"

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Usage: calibrate [-absmem <mem in GB>] [-beam-on-source] [-p <phases.txt> <gains.txt>] [-pf <faraday.txt>] [-px <crossterms.txt>] [-minuv <min uvw dist in m>] [-maxuv <max uvw dist in m>] [-a <min-accuracy> <stop-accuracy>] [-i <niter>] [-j <threads>] [-m <model>] [-scalar] [-diag] [-rhs <rhs solutions>] [-rotation] [-applybeam] [-t timesteps] [-datacolumn <name>] [-quiet] <measurementset.ms> <solutions.bin>\n\n"
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
			savePlotFiles = false, beamOnSource = false, applyBeam = false,
			onlyScalar = false, onlyDiag = false, onlyRotation = false, doQuiet = false;
		std::string plotPhaseFile, plotGainFile, plotFaradayFile, crossTermsPlotFile, modelFile, rhsSolutionFile;
		size_t niter = CalibrationMethod::DefaultNIter(), solutionInterval = 0;
		std::string dataColumnName = "DATA";
		double
			minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
			stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy(),
			minUVW = 0.0,
		  maxUVW = 5000.0,
		  absmem = 0.0;
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
			else if(param == "absmem")
			  {
			    absmem = atof(argv[argi+1]);
			    argi +=2;
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
			else if(param == "maxuv")
			{
				maxUVW = atof(argv[argi+1]);
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
			else if(param == "quiet")
			{
				doQuiet = true;
				++argi;
			}
			else throw std::runtime_error(std::string("Invalid parameter ") + argv[argi]);
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
		calibrator.SetMinUVW(minUVW);
		calibrator.SetMaxUVW(maxUVW);
		calibrator.SetApplyBeam(applyBeam);
		calibrator.SetBeamOnSource(beamOnSource);
		calibrator.SetOnlyScalar(onlyScalar);
		calibrator.SetOnlyDiag(onlyDiag);
		calibrator.SetRHSSolutionFile(rhsSolutionFile);
		calibrator.SetOnlyRotation(onlyRotation);
		calibrator.SetSolutionOutputFilename(outName);
		calibrator.SetVerbose(!doQuiet);
		calibrator.SetSavePlotFiles(savePlotFiles);
		calibrator.SetDataColumnName(dataColumnName);
		calibrator.SetAbsMem(absmem);
		if(savePlotFiles)
		{
			calibrator.SetPlotFilenames(plotPhaseFile, plotGainFile);
		}
		calibrator.Perform();
	}
}
