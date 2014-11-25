#ifndef SOURCE_MATCHER_H
#define SOURCE_MATCHER_H

#include "model.h"
#include "progressbar.h"

class SourceMatcher
{
public:
	enum MatchingType { FindNewMatching, AddSpectraMatching, AvgSpectraMatching };
	
	static void Match(MatchingType matchingType, double distanceInRad, double weight, Model& baseModel, Model& addedModel, Model& outputModel)
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
			size_t ambiguousMatchCount = 0;
			for(size_t i=0; i!=addedModel.SourceCount(); ++i)
			{
				progress.SetProgress(i, addedModel.SourceCount());
				const ModelSource& addSource = addedModel.Source(i);
				
				size_t sourceMatchCount = 0;
				const ModelSource* matchingBaseSource = 0;
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
									++sourceMatchCount;
									matchingBaseSource = &baseSource;
								}
								break;
							}
						}
					}
				}
				
				if(sourceMatchCount > 1)
					ambiguousMatchCount++;
				else if(sourceMatchCount == 1)
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
			
			std::cout << warningsBuffer.str();
			std::cout << "Matched: " << matchedCount << " Ambigious: " << ambiguousMatchCount << " Catalogues: " << baseModel.SourceCount() << " & " << addedModel.SourceCount() << '\n';
		}
	}
};

#endif
