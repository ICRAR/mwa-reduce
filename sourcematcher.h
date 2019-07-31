#ifndef SOURCE_MATCHER_H
#define SOURCE_MATCHER_H

#include "model/model.h"

#include "progressbar.h"
#include "units/angle.h"

#include <fstream>

class SourceMatcher
{
public:
	enum MatchingType { FindNewMatching, AddSpectraMatching, AvgSpectraMatching, CompareMatching };
	
	SourceMatcher() : _acceptAmbigiousDistanceFactor(10), _keepBaseSources(false)
	{ }
	
	void SetAppendAmbigiousFilename(const std::string& filename)
	{
		_appendAmbigiousFilename = filename;
	}
	
	void SetKeepBaseSources(bool keepBaseSources)
	{
		_keepBaseSources = keepBaseSources;
	}
	
	void Match(MatchingType matchingType, double distanceInRad, double weight, Model& baseModel, Model& addedModel, Model& outputModel)
	{
		baseModel.Sort();
		addedModel.Sort();
		
		if(matchingType == FindNewMatching)
		{
			matchNew(distanceInRad, baseModel, addedModel, outputModel);
		}
		else {
			matchCombine(matchingType, distanceInRad, weight, baseModel, addedModel, outputModel);
		}
	}
	
private:
	void matchNew(double distanceInRad, Model& baseModel, Model& addedModel, Model& outputModel)
	{
		size_t matchedCount = 0;
		ProgressBar progress("Matching");
		for(size_t i=0; i!=addedModel.SourceCount(); ++i)
		{
			progress.SetProgress(i, addedModel.SourceCount());
			const ModelSource& addSource = addedModel.Source(i);
			
			bool isMatched = false;
			for(ModelSource::const_iterator compIter=addSource.begin(); compIter!=addSource.end(); ++compIter)
			{
				for(size_t j=0; j!=baseModel.SourceCount(); ++j)
				{
					const ModelSource& baseSource = baseModel.Source(j);
					
					for(ModelSource::const_iterator baseIter=baseSource.begin(); baseIter!=baseSource.end(); ++baseIter)
					{
						double sourceDist =
							ImageCoordinates::AngularDistance(compIter->PosRA(), compIter->PosDec(), baseIter->PosRA(), baseIter->PosDec());
							
						if(sourceDist <= distanceInRad)
						{
							isMatched = true;
							break;
						}
					}
					if(isMatched) break;
				}
				if(isMatched) break;
			}
			
			if(isMatched)
				++matchedCount;
			else {
				outputModel.AddSource(addSource);
			}
		}
		
		progress.SetProgress(1,1);
		
		std::cout << "Unmatched: " << (addedModel.SourceCount()-matchedCount) << " Catalogues: " << baseModel.SourceCount() << " & " << addedModel.SourceCount() << '\n';
	}
	
	void matchCombine(MatchingType matchingType, double distanceInRad, double weight, Model& baseModel, Model& addedModel, Model& outputModel)
	{
		std::unique_ptr<std::ofstream> csvFile;
		if(matchingType == CompareMatching)
		{
			csvFile.reset(new std::ofstream("matched.csv"));
			(*csvFile) << "RA,Dec,S1,S2,S1/S2\n";
		}
		
		std::vector<SourceInfo> infos(baseModel.SourceCount());
		for(size_t i=0; i!=baseModel.SourceCount(); ++i)
			baseModel.Source(i).SetUserData(&infos[i]);
		
		size_t matchedCount = 0;
		ProgressBar progress("Matching");
		std::ostringstream warningsBuffer;
		Model ambigiousModel;
		try {
			if(!_appendAmbigiousFilename.empty())
				ambigiousModel = Model(_appendAmbigiousFilename.c_str());
		} catch(std::exception&) { }
		size_t ambiguousMatchCount = 0;
		for(size_t i=0; i!=addedModel.SourceCount(); ++i)
		{
			const ModelSource& addSource = addedModel.Source(i);
			progress.SetProgress(i, addedModel.SourceCount());
			
			const ModelSource* matchingBaseSource = 0;
			std::vector<const ModelSource*> matchingBaseSources;
			for(ModelSource::const_iterator compIter=addSource.begin(); compIter!=addSource.end(); ++compIter)
			{
				for(size_t j=0; j!=baseModel.SourceCount(); ++j)
				{
					const ModelSource& baseSource = baseModel.Source(j);
					SourceInfo& baseInfo = *static_cast<SourceInfo*>(baseSource.UserData());
					
					for(ModelSource::const_iterator baseIter=baseSource.begin(); baseIter!=baseSource.end(); ++baseIter)
					{
						double sourceDist =
							ImageCoordinates::AngularDistance(compIter->PosRA(), compIter->PosDec(), baseIter->PosRA(), baseIter->PosDec());
							
						if(sourceDist <= distanceInRad && !baseInfo.isMatched)
						{
							if(matchingBaseSource != &baseSource)
							{
								matchingBaseSource = &baseSource;
								matchingBaseSources.push_back(&baseSource);
							}
							break;
						}
					}
				}
			}
			
			if(matchingBaseSources.size() > 1)
			{
				warningsBuffer << addSource.Name() << " ambigious match, with:\n";
				std::vector<std::pair<double, const ModelSource*>> distances;
				for(size_t i=0; i!=matchingBaseSources.size(); ++i)
				{
					const ModelSource& s = *matchingBaseSources[i];
					double sourceDist =
						ImageCoordinates::AngularDistance(s.MeanRA(), s.MeanDec(), addSource.MeanRA(), addSource.MeanDec());
					distances.push_back(std::make_pair(sourceDist, &s));
					warningsBuffer << " - At " << Angle::ToNiceString(sourceDist) << " : " << s.Name() << '\n';
				}
				std::sort(distances.begin(), distances.end());
				if(distances[0].first*_acceptAmbigiousDistanceFactor <= distances[1].first)
				{
					warningsBuffer << "Nearest is " << _acceptAmbigiousDistanceFactor << "x closer, accepting " << distances[0].second->Name() << " as match.\n";
					matchingBaseSource = distances[0].second;
					matchingBaseSources.assign(1, matchingBaseSource);
				}
				else {
					ambiguousMatchCount++;
					ambigiousModel.AddSource(addSource);
				}
			}
			
			if(matchingBaseSources.size() == 1)
			{
				++matchedCount;
				
				SourceInfo& baseInfo = *static_cast<SourceInfo*>(matchingBaseSource->UserData());
				baseInfo.isMatched = true;
				
				ModelSource newSource(*matchingBaseSource);
				if(newSource.ComponentCount() != 1 || addSource.ComponentCount() != 1)
					warningsBuffer << "Warning: matching sources have more than one component, skipping.\n";
				else {
					if(!newSource.front().HasMeasuredSED())
						throw std::runtime_error("Excepting model with measured SEDs");
					if(matchingType == AddSpectraMatching)
					{
						newSource.front().MSED().CombineMeasurements(addSource.front().MSED());
						outputModel.AddSource(newSource);
					}
					else {
						newSource.front().MSED().CombineMeasurementsWithAveraging(addSource.front().MSED(), weight);
						outputModel.AddSource(newSource);
						if(matchingType == CompareMatching)
						{
							long double freq = matchingBaseSource->GetIntegratedMSED().CentreFrequency();
							long double s1 = matchingBaseSource->GetIntegratedMSED().FluxAtFrequencyFromIndex(freq, 0);
							long double s2 = addSource.GetIntegratedMSED().FluxAtFrequencyFromIndex(freq, 0);
							(*csvFile) <<
								newSource.MeanRA()*180.0/M_PI << ',' <<
								newSource.MeanDec()*180.0/M_PI << ',' <<
								s1 << ',' <<
								s2 << ',' <<
								(s1/s2) << '\n';
						}
					}
				}
			}
		}
		
		if(_keepBaseSources)
		{
			for(Model::const_iterator s=baseModel.begin();
					s!=baseModel.end(); ++s)
			{
				SourceInfo& baseInfo = *static_cast<SourceInfo*>(s->UserData());
				if(!baseInfo.isMatched)
				{
					outputModel.AddSource(*s);
				}
			}
		}
		
		progress.SetProgress(1,1);
		
		if(!_appendAmbigiousFilename.empty())
			ambigiousModel.Save(_appendAmbigiousFilename);
		
		std::cout << warningsBuffer.str();
		std::cout << "Matched: " << matchedCount << " Ambigious: " << ambiguousMatchCount << " Catalogues: " << baseModel.SourceCount() << " & " << addedModel.SourceCount() << '\n';
	}

	struct SourceInfo {
		SourceInfo() : isMatched(false) { }
		bool isMatched;
	};
	double _acceptAmbigiousDistanceFactor;
	std::string _appendAmbigiousFilename;
	bool _keepBaseSources;
};

#endif
