#include "uvector.h"
#include "model.h"
#include "ionsolutionfile.h"
#include "beamevaluator.h"
#include "banddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <iostream>

/**
For each measurement set:

- Calculate apparent flux for source
- Multiple by peeled gain
- Weight by beam

- Add source flux*weight to source sum, add weight to source weight

Finally:
- Divide by weight^2
- Convert to Stokes
*/

std::map<std::string, size_t> sourceNameToIndex;
ao::uvector<std::complex<double>> gainPerSource, weightsPerSource;
double centralFrequency = 0.0;

void addSolutionFile(Model& model, const char* solutionFilename, const char* msFilename)
{
	std::cout << "Opening " << solutionFilename << "... " << std::flush;
	IonSolutionFile file;
	file.OpenForReading(solutionFilename);
	
	std::vector<std::vector<size_t>> sourceIndexPerCluster;
	for(size_t direction=0; direction!=file.DirectionCount(); ++direction)
	{
		std::string clusterName;
		std::vector<std::string> sourceNames;
		file.ReadClusterMetaInfo(clusterName, sourceNames);
		if(file.PolarizationCount() != 1)
			throw std::runtime_error("This solution file format is unknown");
		
		std::vector<size_t> sourceIndices(sourceNames.size());
		for(size_t s=0; s!=sourceNames.size(); ++s)
		{
			std::map<std::string, size_t>::const_iterator nameToIndexIter = sourceNameToIndex.find(sourceNames[s]);
			if(nameToIndexIter == sourceNameToIndex.end())
				throw std::runtime_error("Solutions contained source which is not in the model");
			sourceIndices[s] = nameToIndexIter->second;
		}
		sourceIndexPerCluster.push_back(sourceIndices);
	}
	
	casa::MeasurementSet ms(msFilename);
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	casa::MEpoch
		firstTime = timeColumn(0),
		lastTime = timeColumn(ms.nrow()-1),
		time = casa::MEpoch(casa::MVEpoch(0.5 * (firstTime.getValue().get() + lastTime.getValue().get())), firstTime.getRef());
	BeamEvaluator beamEvaluator(ms, false);
	beamEvaluator.SetTime(time);
	BandData bandData(ms.spectralWindow());
	if(centralFrequency == 0.0)
		centralFrequency = bandData.CentreFrequency();
	else {
		if(bandData.CentreFrequency() != centralFrequency)
			throw std::runtime_error("Measurement sets have different frequency coverage");
	}
	
	char hasBeenProcessedMark = 1, hasNotBeenProcessedMark = 0;
	for(Model::iterator s=model.begin(); s!=model.end(); ++s)
		s->SetUserData(&hasNotBeenProcessedMark);
	
	double gainSum = 0.0;
	for(size_t direction=0; direction!=file.DirectionCount(); ++direction)
	{
		double gain = file.ReadAverageSolution(IonSolutionFile::GainSolution, 0, direction);
		
		if(std::isfinite(gain))
		{
			gainSum += gain;
			
			const std::vector<size_t>& sourceIndices = sourceIndexPerCluster[direction];
			
			for(size_t s=0; s!=sourceIndices.size(); ++s)
			{
				size_t index = sourceIndices[s];
				ModelSource& source = model.Source(index);
				source.SetUserData(&hasBeenProcessedMark);
				std::complex<double> temp[4], beamSq[4];
				beamEvaluator.EvaluateAbsToApparentGain(source.MeanRA(), source.MeanDec(), centralFrequency, temp);
				
				// (w is assumed unity for now)
				// gainsum += B* w g B
				Matrix2x2::HermATimesB(beamSq, temp, temp);
				Matrix2x2::MultiplyAdd(&gainPerSource[index*4], beamSq, gain);
				
				// weightsum += (B* B) w
				//Matrix2x2::ATimesB(temp, beamSq, beamSq);
				Matrix2x2::Add(&weightsPerSource[index*4], beamSq);
			}
		}
	}
	
	for(size_t index=0; index!=model.SourceCount(); ++index)
	{
		ModelSource& source = model.Source(index);
		if(source.UserData() == &hasNotBeenProcessedMark)
		{
			std::complex<double> temp[4], beamSq[4];
			beamEvaluator.EvaluateAbsToApparentGain(source.MeanRA(), source.MeanDec(), centralFrequency, temp);
			// Only add to weight, since it hasn't been peeled
			Matrix2x2::HermATimesB(beamSq, temp, temp);
			Matrix2x2::Add(&weightsPerSource[index*4], beamSq);
		}
	}

	std::cout << "Avg gain=" << (gainSum / file.DirectionCount()) << '\n';
}

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout << "Syntax: collectionsolutions <full model> <output model> [{<solution file1> <ms>} ...]\n";
		return 0;
	}
	int argi = 1;
	const char* modelFilename = argv[argi];
	const char* outputModelFilename = argv[argi+1];
	
	std::cout << "Reading model...\n";
	Model model(modelFilename);
	
	gainPerSource.assign(model.SourceCount()*4, 0.0);
	weightsPerSource.assign(model.SourceCount()*4, 0.0);
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const std::string& name = model.Source(i).Name();
		if(sourceNameToIndex.count(name) != 0)
			throw std::runtime_error(std::string("Double source in model: ") + name);
		sourceNameToIndex.insert(std::make_pair(name, i));
	}
	
	argi+= 2;
	while(argi < argc)
	{
		addSolutionFile(model, argv[argi], argv[argi+1]);
		argi += 2;
	}
	
	Model outputModel;
	for(size_t s=0; s!=model.SourceCount(); ++s)
	{
		std::complex<double> *sw = &weightsPerSource[s*4];
		if(sw[0] != 0.0 || sw[1] != 0.0 || sw[2] != 0.0 || sw[3] != 0.0)
		{
			std::complex<double> *sg = &gainPerSource[s*4];
			std::complex<double> correctedLinear[4];
			
			// Calculate: W*GW sum(W*W)^-1
			// sum(W*W) is in variable sw, and W*GW in sg.
			if(Matrix2x2::Invert(sw))
			{
				Matrix2x2::ATimesB(correctedLinear, sg, sw);
				double avgGain = (correctedLinear[0].real() + correctedLinear[3].real()) * 0.5;
				ModelSource newSource(model.Source(s));
				newSource *= avgGain;
				std::cout << newSource.Name() << " " << avgGain << '\n';
				outputModel.AddSource(newSource);
			}
		}
	}
	outputModel.Save(outputModelFilename);
}
