#ifndef SOURCE_MATCHER_H
#define SOURCE_MATCHER_H

#include "model.h"
#include "progressbar.h"
#include "angle.h"

class SourceMatcher
{
public:
	enum MatchingType { FindNewMatching, AddSpectraMatching, AvgSpectraMatching };
	
	SourceMatcher() : _acceptAmbigiousDistanceFactor(10)
	{ }
	
	void SetAppendAmbigiousFilename(const std::string& filename)
	{
		_appendAmbigiousFilename = filename;
	}
	
	void Match(MatchingType matchingType, double distanceInRad, double weight, Model& baseModel, Model& addedModel, Model& outputModel)
	{
		baseModel.SortOnBrightness();
		addedModel.SortOnBrightness();
		
		size_t matchedCount = 0;
		
		std::ostringstream warningsBuffer;
		ProgressBar progress("Matching");
		
		if(matchingType == FindNewMatching)
		{
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
		else {
			Model ambigiousModel;
			try {
				if(!_appendAmbigiousFilename.empty())
					ambigiousModel = Model(_appendAmbigiousFilename.c_str());
			} catch(std::exception&) { }
			size_t ambiguousMatchCount = 0;
			for(size_t i=0; i!=addedModel.SourceCount(); ++i)
			{
				progress.SetProgress(i, addedModel.SourceCount());
				const ModelSource& addSource = addedModel.Source(i);
				
				const ModelSource* matchingBaseSource = 0;
				std::vector<const ModelSource*> matchingBaseSources;
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
					ModelSource newSource(*matchingBaseSource);
					if(newSource.ComponentCount() != 1 || addSource.ComponentCount() != 1)
						warningsBuffer << "Warning: matching sources have more than one component, skipping.\n";
					else {
						if(matchingType == AddSpectraMatching)
							newSource.front().SED().CombineMeasurements(addSource.front().SED());
						else
							newSource.front().SED().CombineMeasurementsWithAveraging(addSource.front().SED(), weight);
						outputModel.AddSource(newSource);
					}
				}
			}
			progress.SetProgress(1,1);
			
			if(!_appendAmbigiousFilename.empty())
				ambigiousModel.Save(_appendAmbigiousFilename);
			
			std::cout << warningsBuffer.str();
			std::cout << "Matched: " << matchedCount << " Ambigious: " << ambiguousMatchCount << " Catalogues: " << baseModel.SourceCount() << " & " << addedModel.SourceCount() << '\n';
		}
	}
	
private:
	double _acceptAmbigiousDistanceFactor;
	std::string _appendAmbigiousFilename;
};

#endif
