#include "model.h"
#include "progressbar.h"

int main(int argc, char* argv[])
{
	if(argc < 6)
	{
		std::cout <<
			"Syntax: matchsources <type> [options] <base model> <additional model> <arcmin distance> <output model>\n"
			"types:\n"
			" - findnew    -- Find any sources in additional and not in base\n"
			" - addspectra -- Match sources and combine frequency info of both catalogues\n"
			" - avgspectra -- Match sources and combine frequency info of both catalogues. Duplicate measurements are averaged.\n";
			return 0;
	}
	
	enum { FindNewMatching, AddSpectraMatching, AvgSpectraMatching } matchingType;
	
	const std::string matchingTypeStr(argv[1]);
	
	if(matchingTypeStr == "findnew")
		matchingType = FindNewMatching;
	else if(matchingTypeStr == "addspectra")
		matchingType = AddSpectraMatching;
	else if(matchingTypeStr == "avgspectra")
		matchingType = AvgSpectraMatching;
	else
		throw std::runtime_error("Unknown matching type specified");
	
	double weight = 0.5;
	size_t argi = 2;
	while(argv[argi][0] == '-')
	{
		std::string p(&argv[argi][1]);
		if(p == "weight")
		{
			weight = atof(argv[argi+1]);
			++argi;
		}
		else throw std::runtime_error("Invalid parameter");
		++argi;
	}
	
	Model baseModel(argv[argi]), addModel(argv[argi+1]), restModel;
	double distance = atof(argv[argi+2])*(M_PI/60.0/180.0);
	const std::string restModelFilename(argv[argi+3]);
	
	baseModel.SortOnBrightness();
	addModel.SortOnBrightness();
	
	size_t matchedCount = 0;
	
	ProgressBar progress("Matching");
	
	if(matchingType == FindNewMatching)
	{
		for(size_t i=0; i!=addModel.SourceCount(); ++i)
		{
			progress.SetProgress(i, addModel.SourceCount());
			const ModelSource& addSource = addModel.Source(i);
			
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
							
						if(sourceDist <= distance)
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
				restModel.AddSource(addSource);
			}
		}
		
		progress.SetProgress(1,1);
		
		std::cout << "Unmatched: " << (addModel.SourceCount()-matchedCount) << " Catalogues: " << baseModel.SourceCount() << " & " << addModel.SourceCount() << '\n';
	}
	else {
		size_t ambiguousMatchCount = 0;
		for(size_t i=0; i!=addModel.SourceCount(); ++i)
		{
			progress.SetProgress(i, addModel.SourceCount());
			const ModelSource& addSource = addModel.Source(i);
			
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
							
						if(sourceDist <= distance)
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
					std::cout << "Warning: matching sources have more than one component, skipping.\n";
				else {
					if(matchingType == AddSpectraMatching)
						newSource.front().SED().CombineMeasurements(addSource.front().SED());
					else
						newSource.front().SED().CombineMeasurementsWithAveraging(addSource.front().SED(), weight);
					restModel.AddSource(newSource);
				}
			}
		}
		progress.SetProgress(1,1);
		
		std::cout << "Matched: " << matchedCount << " Ambigious: " << ambiguousMatchCount << " Catalogues: " << baseModel.SourceCount() << " & " << addModel.SourceCount() << '\n';
	}
	
	
	restModel.Save(restModelFilename);
}
