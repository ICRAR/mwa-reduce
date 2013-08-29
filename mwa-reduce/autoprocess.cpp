#include <iostream>

#include "beamevaluator.h"
#include "model.h"
#include "banddata.h"
#include "imagecoordinates.h"
#include "peeler.h"

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
		std::cout << "Syntax: autoprocess [options] <catalogue> <ms>\n";
		return 1;
	}
	int argi = 1;
	bool doExecute = false;
	while(argv[argi][0] == '-')
	{
		std::string param(&argv[argi][1]);
		if(param == "go")
		{
			doExecute = true;
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
	
	BeamEvaluator beamEvaluator(ms, false);
	
	std::vector<std::pair<double, ModelSource*>> sources;
	
	for(Model::iterator srcIter=catalogue.begin(); srcIter!=catalogue.end(); ++srcIter)
	{
		ModelSource& source = *srcIter;
		std::complex<double> fluxMatrix[4];
		for(size_t p=0; p!=4; ++p)
			fluxMatrix[p] = source.SED().IntegratedFlux(bandData.LowestFrequency(), bandData.HighestFrequency(), p);
		beamEvaluator.AbsToApparent(source.PosRA(), source.PosDec(), fluxMatrix);
		double fluxStokesI = (fluxMatrix[0].real() + fluxMatrix[3].real()) * 0.5;
		sources.push_back(std::make_pair(fluxStokesI, &source));
	}
	
	std::sort(sources.rbegin(), sources.rend());

	std::vector<ModelSource*> calibrateSources, peelSources, subtractSources;
	
	std::cout << "Strongest apparent sources:\n";
	for(size_t i=0; i!=sources.size(); ++i)
	{
		std::pair<double, ModelSource*>& src = sources[i];
		double distanceDeg = (180.0 / M_PI) *
			ImageCoordinates::AngularDistance(phaseCentreRA, phaseCentreDec, src.second->PosRA(), src.second->PosDec());
		double
			distanceNice = round(distanceDeg*10.0)*0.1,
			fluxNice = round(src.first*10.0)*0.1;
		std::cout << src.second->Name() << " (" << fluxNice << " Jy/beam, distance=" << distanceNice << " deg)\n";
		
		// Determine what to do with it
		if(src.first >= peelThreshold)
		{
			if(distanceDeg <= maxCalibrateDist)
			{
				calibrateSources.push_back(src.second);
			}
			else {
				double runnerUpFlux = 0.0;
				if(i+1 < sources.size())
					runnerUpFlux = sources[i+1].first;
				if(runnerUpFlux * peelMinRunnerUpFactor > src.first)
				{
					std::cout << " -> Runner-up source " << sources[i+1].second->Name() << " is too bright for automated peeling of this source.\n";
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
		std::cout << "- No self-calibration (no obvious strong sources close to the centre).\n";
	else
		std::cout << "- Self-calibrate using " << sourceList(calibrateSources) << '\n';
	
	if(peelSources.empty())
		std::cout << "- No peeling.\n";
	else
		std::cout << "- Peel out " << sourceList(peelSources) << '\n';
	
	if(subtractSources.empty())
		std::cout << "- No spectral source subtraction.\n";
	else
		std::cout << "- Spectrally subtract " << sourceList(subtractSources) << '\n';
	
	std::cout << '\n';	
	if(calibrateSources.size() + peelSources.size() + subtractSources.size() == 0)
		std::cout << "No tasks to perform.\n";
	else if(!doExecute)
		std::cout << "Not performing tasks: specify '-go' on command line to execute.\n";
	else
	{
		if(!peelSources.empty())
		{
			Model peelModel;
			for(std::vector<ModelSource*>::const_iterator i=peelSources.begin(); i!=peelSources.end(); ++i)
				peelModel.AddSource(**i);
			
			Peeler peeler(ms);
			
			peeler.SetNIter(1000);
			peeler.SetLimit(0.001);
			peeler.SetModel(peelModel);
			peeler.SetDataColumName("CORRECTED_DATA");
			peeler.SetSolutionInterval(4);
			
			peeler.Perform();
		}
	}
}
