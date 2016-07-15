#include "model/model.h"

#include "sourcematcher.h"

int main(int argc, char* argv[])
{
	if(argc < 6)
	{
		std::cout <<
			"Syntax: matchsources <type> [options] <output model> <arcmin distance> <base model> <additional model> [models...]\n"
			"types:\n"
			" - findnew    -- Find any sources in additional and not in base\n"
			" - addspectra -- Match sources and combine frequency info of both catalogues\n"
			" - avgspectra -- Match sources and combine frequency info of both catalogues. Duplicate measurements are averaged.\n";
			return 0;
	}
	
	SourceMatcher::MatchingType matchingType;
	const std::string matchingTypeStr(argv[1]);
	
	if(matchingTypeStr == "findnew")
		matchingType = SourceMatcher::FindNewMatching;
	else if(matchingTypeStr == "addspectra")
		matchingType = SourceMatcher::AddSpectraMatching;
	else if(matchingTypeStr == "avgspectra")
		matchingType = SourceMatcher::AvgSpectraMatching;
	else if(matchingTypeStr == "compare")
		matchingType = SourceMatcher::CompareMatching;
	else
		throw std::runtime_error("Unknown matching type specified");
	
	double weight = 0.5;
	size_t argi = 2;
	SourceMatcher matcher;
	while(argv[argi][0] == '-')
	{
		std::string p(&argv[argi][1]);
		if(p == "weight")
		{
			++argi;
			weight = atof(argv[argi]);
		}
		else if(p == "append-ambigious")
		{
			++argi;
			matcher.SetAppendAmbigiousFilename(argv[argi]);
		}
		else if(p == "keep-base")
		{
			matcher.SetKeepBaseSources(true);
		}
		else throw std::runtime_error("Invalid parameter");
		++argi;
	}
	
	double distanceInRad = atof(argv[argi+1])*(M_PI/60.0/180.0);
	const std::string restModelFilename(argv[argi]);
	std::cout << "Reading " << argv[argi+2] << "...\n";
	Model baseModel(argv[argi+2]);
	argi+=3;
	do {
		std::string addedModelFilename(argv[argi]);
		Model restModel;
		std::cout << "Reading " << addedModelFilename << "...\n";
		Model addedModel(addedModelFilename);
		matcher.Match(matchingType, distanceInRad, weight, baseModel, addedModel, restModel);
		
		baseModel = restModel;
		++argi;
	} while(argi < size_t(argc));
	
	baseModel.Save(restModelFilename);
}
