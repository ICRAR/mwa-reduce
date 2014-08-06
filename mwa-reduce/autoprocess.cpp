#include <iostream>

#include "beamevaluator.h"
#include "model.h"
#include "banddata.h"
#include "imagecoordinates.h"
#include "peeler.h"
#include "calibrator.h"
#include "solutionapplier.h"
#include "subtractor.h"
#include "spectrumsubtractor.h"
#include "calibrationmethod.h"

std::string sourceList(const std::vector<ModelSource*>& sources)
{
	std::ostringstream str;
	if(sources.size() > 1)
		str << "sources ";
	else
		str << "source ";
	
	str << sources[0]->Name();
	for(size_t i=1; i != sources.size(); ++i)
		str << ", " << sources[i]->Name();
	
	return str.str();
}

int main(int argc, char* argv[])
{
	if(argc < 2) {
		std::cout << "Syntax:\n\t autoprocess [options] <catalogue> <ms>\n\n"
			"Options:\n"
			"-go\n"
			"  Execute: without it, nothing will be done, only suggestions will be given.\n"
			"-noselfcal\n"
			"  Disable self cal.\n"
			"-subtract\n"
			"  Enable direct subtraction of sources.\n"
			"-nosubtract\n"
			"  Disable direct subtraction of sources.\n"
			"-nopeel\n"
			"   Disable peeling of sources.\n"
			"-peel-with-runnerups\n"
			"   Enable peeling even when runner ups are strong.\n"
			"-niter <iterations>\n"
			"   Number of selfcal/peeling iterations.\n"
			"-save-solution-files\n"
			"   Save the solution files, such as 'peel-sol-ant25.txt', etc.\n"
			"-a <min accuracy> <stop accuracy>\n"
			"   Set accuracy at which to accept or immediately accept a solution.\n";
		return 1;
	}
	int argi = 1;
	bool doExecute = false, doSelfCal = false, doPeel = true, doSubtract = false, verboseOnPolarizations = false,
		noPeelingRunnerupCheck = false, saveSolutionFiles = false;
	double
		minAccuracy = CalibrationMethod::DefaultMinAccuracy(),
		stopAccuracy = CalibrationMethod::DefaultStoppingAccuracy();
	size_t nIter = CalibrationMethod::DefaultNIter();
	while(argv[argi][0] == '-')
	{
		std::string param(&argv[argi][1]);
		if(param == "go")
		{
			doExecute = true;
		}
		else if(param == "noselfcal")
		{
			doSelfCal = false;
		}
		else if(param == "nopeel")
		{
			doPeel = false;
		}
		else if(param == "subtract")
		{
			doSubtract = true;
		}
		else if(param == "nosubtract")
		{
			doSubtract = false;
		}
		else if(param == "peel-with-runnerups")
		{
			noPeelingRunnerupCheck = true;
		}
		else if(param == "niter")
		{
			++argi;
			nIter = atoi(argv[argi]);
		}
		else if(param == "save-solution-files")
		{
			saveSolutionFiles = true;
		}
		else if(param == "vp")
		{
			verboseOnPolarizations = true;
		}
		else if(param == "a")
		{
			minAccuracy = atof(argv[argi+1]);
			stopAccuracy = atof(argv[argi+2]);
			argi+=2;
		}
		else throw std::runtime_error("Invalid parameter");
		++argi;
	}
	
	Model catalogue(argv[argi]);
	casa::MeasurementSet ms(argv[argi+1]);
	BandData bandData(ms.spectralWindow());
	
	casa::MSField fieldTable = ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	
	double
		subtractThreshold = 20.0,
		peelThreshold = 50.0,
		peelMinRunnerUpFactor = 3.0,
		maxCalibrateDist = 12.5;
	
	bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
	std::string dataColumn;
	if(hasCorrected) {
		std::cout << "This measurement set has corrected data: tasks will be applied on the corrected data column.\n";
		dataColumn = "CORRECTED_DATA";
	} else {
		std::cout << "No corrected data in set: tasks will be applied on the data column.\n";
		dataColumn= "DATA";
	}
		
	BeamEvaluator beamEvaluator(ms, false);
	
	std::vector<std::pair<double, ModelSource*>> sources;
	
	for(Model::iterator srcIter=catalogue.begin(); srcIter!=catalogue.end(); ++srcIter)
	{
		ModelSource& source = *srcIter;
		std::complex<double> fluxLinMatrix[4];
		double fluxStokesMatrix[4];
		for(size_t p=0; p!=4; ++p)
			fluxStokesMatrix[p] = source.TotalFlux(bandData.LowestFrequency(), bandData.HighestFrequency(), Polarization::IndexToStokes(p));
		Polarization::StokesToLinear(fluxStokesMatrix, fluxLinMatrix);
		beamEvaluator.AbsToApparent(source.Peak().PosRA(), source.Peak().PosDec(), fluxLinMatrix);
		double fluxStokesI = (fluxLinMatrix[0].real() + fluxLinMatrix[3].real()) * 0.5;
		sources.push_back(std::make_pair(fluxStokesI, &source));
	}
	
	std::sort(sources.rbegin(), sources.rend());

	std::vector<ModelSource*> calibrateSources, peelSources, subtractSources;
	bool peelingSourceSkipped = false;
	
	std::cout << "Strongest apparent sources:\n";
	for(size_t i=0; i!=sources.size(); ++i)
	{
		std::pair<double, ModelSource*>& src = sources[i];
		double distanceDeg = (180.0 / M_PI) *
			ImageCoordinates::AngularDistance(phaseCentreRA, phaseCentreDec, src.second->Peak().PosRA(), src.second->Peak().PosDec());
		double
			distanceNice = round(distanceDeg*10.0)*0.1,
			fluxNice = round(src.first*10.0)*0.1;
		
		std::cout << src.second->Name() << " (" << fluxNice << " Jy/beam, distance=" << distanceNice << " deg)\n";
		if(verboseOnPolarizations) {
			std::complex<double> fluxLinMatrix[4];
			double fluxStokesMatrix[4];
			const ModelSource& source = *src.second;
			for(size_t p=0; p!=4; ++p)
				fluxStokesMatrix[p] = source.TotalFlux(bandData.LowestFrequency(), bandData.HighestFrequency(), Polarization::IndexToStokes(p));
			Polarization::StokesToLinear(fluxStokesMatrix, fluxLinMatrix);
			beamEvaluator.AbsToApparent(source.Peak().PosRA(), source.Peak().PosDec(), fluxLinMatrix);
			std::cout << round(fluxLinMatrix[0].real()*10.0)*0.1;
			for(size_t p=1; p!=4; ++p)
				std::cout << ' ' << round(fluxLinMatrix[p].real()*10.0)*0.1;
			std::cout << '\n';
		}
		
		// Determine what to do with it
		if(src.first >= peelThreshold && !peelingSourceSkipped)
		{
			if(distanceDeg <= maxCalibrateDist && peelSources.empty() && doSelfCal)
			{
				calibrateSources.push_back(src.second);
			}
			else {
				double runnerUpFlux = 0.0;
				if(i+1 < sources.size())
					runnerUpFlux = sources[i+1].first;
				if(runnerUpFlux * peelMinRunnerUpFactor > src.first && !noPeelingRunnerupCheck)
				{
					std::cout << " -> Runner-up source " << sources[i+1].second->Name() << " is too bright for automated peeling of this source.\n";
					peelingSourceSkipped = true;
				}
				else {
					peelSources.push_back(src.second);
				}
			}
		}
		else if(src.first >= subtractThreshold)
		{
			subtractSources.push_back(src.second);
		}
	}

	std::cout << "\nTasks:\n";
	
	if(calibrateSources.empty())
		std::cout << "- No self-calibration (no strong sources close to the centre).\n";
	else
		std::cout << "- Self-calibrate using " << sourceList(calibrateSources) << " and subtract " << ((calibrateSources.size()>1) ? "these sources.\n" : "this source.\n");
	
	if(peelSources.empty())
		std::cout << "- No peeling.\n";
	else {
		if(doPeel)
			std::cout << "- Peel out " << sourceList(peelSources) << '\n';
		else
			std::cout << "- Advice is to peel out " << sourceList(peelSources) << ", but peeling is disabled.\n";
	}
	
	if(subtractSources.empty())
		std::cout << "- No spectral source subtraction.\n";
	else {
		if(doSubtract)
			std::cout << "- Spectrally subtract " << sourceList(subtractSources) << '\n';
		else
			std::cout << "- Advice is to subtract " << sourceList(subtractSources) << ", but subtraction is disabled.\n";
	}
	
	std::cout << '\n';	
	if(calibrateSources.size() + peelSources.size() + subtractSources.size() == 0)
		std::cout << "No tasks to perform.\n";
	else if(!doExecute)
		std::cout << "Not performing tasks: specify '-go' on command line to execute.\n";
	else
	{
		Model restorationModel;
		
		if(!calibrateSources.empty())
		{
			std::cout << "Calibrating using " << sourceList(calibrateSources) << "...\n";
			
			Calibrator calibrator(ms);
			
			Model calModel;
			for(std::vector<ModelSource*>::const_iterator i=calibrateSources.begin(); i!=calibrateSources.end(); ++i)
			{
				calModel.AddSource(**i);
				restorationModel.AddSource(**i);
			}
			calibrator.SetModel(calModel);
			calibrator.SetDataColumnName(dataColumn);
			calibrator.SetSolutionInterval(0);
			calibrator.SetAccuracy(minAccuracy, stopAccuracy);
			calibrator.SetApplyBeam(true);
			calibrator.SetNIter(nIter);
			calibrator.SetVerbose(true);
			
			calibrator.Perform();
			
			std::cout << "Applying solutions...\n";
			SolutionApplier applier;
			applier.Apply(ms, calibrator.GetSolutionFile());
			
			std::cout << "Subtracting model...\n";
			Subtractor subtractor;
			subtractor.SetApplyBeam(true);
			subtractor.Subtract(ms, calModel);
		}
		
		if(!peelSources.empty() && doPeel)
		{
			std::cout << "Peeling " << sourceList(peelSources) << "...\n";
			
			for(std::vector<ModelSource*>::const_iterator i=peelSources.begin(); i!=peelSources.end(); ++i)
			{
				Model peelModel;
				restorationModel.AddSource(**i);
				ModelSource peelSource = **i;
				std::cout << "- Source " << peelSource.Name() << "\n";
				// Correct for the beam; this is not necessarily as gains are fitted, but will (a)
				// give better initial conditions and (b) the reported gains are true differential gains.
				for(ModelSource::iterator i=peelSource.begin(); i!=peelSource.end(); ++i)
				{
					std::complex<double> beamMatrix[4], beamGain[4];
					beamEvaluator.EvaluateAbsToApparentGain(i->PosRA(), i->PosDec(), beamMatrix);
					Matrix2x2::ATimesHermB(beamGain, beamMatrix, beamMatrix);
					double gain = (beamGain[0].real() + beamGain[3].real()) * 0.5;
					i->SED() *= gain;
				}
				peelModel.AddSource(peelSource);
				
				Peeler peeler(ms);
				
				peeler.SetModel(peelModel);
				peeler.SetDataColumnName(dataColumn);
				peeler.SetSolutionInterval(4);
				peeler.SetAccuracy(minAccuracy, stopAccuracy);
				peeler.SetNIter(nIter);
				peeler.SetSaveSolutionFiles(saveSolutionFiles);
				
				peeler.Perform();
			}
		}
		
		if(!subtractSources.empty() && doSubtract)
		{
			std::cout << "Spectrally subtracting " << sourceList(subtractSources) << "...\n";
			
			Model subtractModel;
			for(std::vector<ModelSource*>::const_iterator src=subtractSources.begin(); src!=subtractSources.end(); ++src)
			{
				ModelSource subtractSource = **src;
				for(ModelSource::iterator i=subtractSource.begin(); i!=subtractSource.end(); ++i)
				{
					std::complex<double> beamMatrix[4], beamGain[4];
					beamEvaluator.EvaluateAbsToApparentGain(i->PosRA(), i->PosDec(), beamMatrix);
					Matrix2x2::ATimesHermB(beamGain, beamMatrix, beamMatrix);
					double gain = (beamGain[0].real() + beamGain[3].real()) * 0.5;
					i->SED() *= gain;
				}
				subtractModel.AddSource(subtractSource);
			}
			
			SpectrumSubtractor subtractor(ms, subtractModel);
			subtractor.SetFittingInterval(4);
			subtractor.SetDataColumnName(dataColumn);
			subtractor.Perform();
			
			// Add the fitted sources to the restoration model
			const Model& fittedModel = subtractor.RestorationModel();
			for(Model::const_iterator src=fittedModel.begin(); src!=fittedModel.end(); ++src)
			{
				restorationModel.AddSource(*src);
			}
		}
		
		restorationModel.Save("model-restore.txt");
		std::cout << "Restoration model containing " << restorationModel.SourceCount() << " sources written to model-restore.txt .\n";
	}
}
