#include "sourcematcher.h"
#include "sedanalyser.h"

#include <iostream>
#include <fstream>
#include <random>

void printRMS(const Model& model)
{
	std::cout << "Running analyser...\n";
	SEDAnalyser sedAnalyser;
	sedAnalyser.Set(model);
	sedAnalyser.Process();
	std::cout << "RMS: " << sedAnalyser.BestRMS() << '\n';
}

struct SEDSet
{
	std::string bandname, shortname, filename;
	double rms, integrationTime;
	Model model;
	bool operator<(const SEDSet& other) const { return filename < other.filename; }
};

void processList(std::vector<SEDSet> sets, Model& fullModel, std::ofstream& logFile)
{
	bool useRMS = true;
	std::set<std::string> bandNames;
	for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
		bandNames.insert(i->bandname);
				
	fullModel = Model();
	double distance = 0.3 * (M_PI/60.0/180.0);
	
	std::vector<Model> bandModels(bandNames.size());
	std::vector<Model>::iterator bandModelIter = bandModels.begin();
	for(std::set<std::string>::const_iterator bandIter=bandNames.begin(); bandIter!=bandNames.end(); ++bandIter,++bandModelIter)
	{
		Model& bandModel = *bandModelIter;
		std::vector<SEDSet> bandSets;
		for(std::vector<SEDSet>::const_iterator setIter=sets.begin(); setIter!=sets.end(); ++setIter)
		{
			if(setIter->bandname == *bandIter)
				bandSets.push_back(*setIter);
		}
		std::cout << "Processing " << bandSets.size() << " sets in band " << *bandIter << '\n';
		
		std::vector<SEDSet>::const_iterator bandSetIter = bandSets.begin();
		std::cout << "Loading " << bandSetIter->shortname << "...\n";
		bandModel = Model(bandSetIter->filename.c_str());
		double weightSum = useRMS ? 1.0/(bandSetIter->rms*bandSetIter->rms) : 1.0;
		
		++bandSetIter;
		while(bandSetIter != bandSets.end())
		{
			double thisWeight = useRMS ? 1.0/(bandSetIter->rms*bandSetIter->rms) : 1.0;
			
			double addInWeight = thisWeight / (weightSum + thisWeight);
			weightSum += thisWeight;
			
			std::cout << "Loading " << bandSetIter->shortname << "...\n";
			Model addedModel(bandSetIter->filename.c_str()), outputModel;
			
			std::cout << "Averaging in " << bandSetIter->shortname << "... (" << bandSetIter->rms << " mJy RMS, " << round(addInWeight*1000.0)/10.0 << "%)\n";
			SourceMatcher matcher;
			matcher.Match(SourceMatcher::AvgSpectraMatching, distance, addInWeight, bandModel, addedModel, outputModel);
			bandModel = outputModel;
			++bandSetIter;
			
			//std::cout << "Band ";
			//printRMS(bandModel);
		}
		
		std::cout << "Theoretical noise in band: " << sqrt(1.0/weightSum) << " mJy.\n";
		bandModel.Save(std::string(*bandIter) + "-band-model.txt");
		if(bandIter == bandNames.begin())
			fullModel = bandModel;
		else {
			std::cout << "Combining bands...\n";
			Model outputModel;
			SourceMatcher matcher;
			matcher.Match(SourceMatcher::AvgSpectraMatching, distance, 0.5, fullModel, bandModel, outputModel);
			fullModel = outputModel;
		}
	}
	
	std::cout << "Running analyser...\n";
	SEDAnalyser sedAnalyser;
	sedAnalyser.Set(fullModel);
	sedAnalyser.Process();
	
	std::vector<SEDSet>::const_iterator i=sets.begin();
	logFile << sedAnalyser.BestRMS() << '\t' << i->shortname;
	++i;
	for(; i!=sets.end(); ++i)
		logFile << ',' << i->shortname;
	logFile << std::endl;
}

void makePermutation(size_t setPerBand, size_t maxSetsInBand, const std::vector<SEDSet>& sets, std::vector<SEDSet>& permutedSets)
{
	permutedSets.clear();
	std::set<std::string> bandNames;
	for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
		bandNames.insert(i->bandname);
	
	std::random_device rd;
	std::mt19937 rng(rd());
	
	for(std::set<std::string>::const_iterator band=bandNames.begin(); band!=bandNames.end(); ++band)
	{
		size_t setCount = 0;
		for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
			if(i->bandname == *band) ++setCount;
			
		std::set<size_t> selectedIndices;
		std::uniform_int_distribution<uint32_t> uint_dist(0, setCount-1);
		
		while(selectedIndices.size() < std::min(setCount, setPerBand))
		{
			size_t index;
			do {
				index = uint_dist(rng);
			} while(selectedIndices.count(index) != 0);
			selectedIndices.insert(index);
			size_t browse = 0;
			for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
			{
				if(i->bandname == *band)
				{
					if(browse == index)
						permutedSets.push_back(*i);
					++browse;
				}
			}
		}
	}
	std::sort(permutedSets.begin(), permutedSets.end());
}

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout <<
			"Syntax: sedcombine <listfile> <outputmodel>\n\n"
			"format of rows in listfile:\n"
			"bandname RMS integrationTime shortname filename\n";
	}
	else {
		std::ifstream listfile(argv[1]);
		std::string outputModelName(argv[2]);
		std::vector<SEDSet> sets;
		
		while(listfile.good())
		{
			SEDSet sedSet;
			listfile >> sedSet.bandname >> sedSet.rms >> sedSet.integrationTime >> sedSet.shortname >> sedSet.filename;
			if(listfile.good())
			{
				sets.push_back(sedSet);
			}
		}
		
		// Count nr of sets in the bands
		std::map<std::string, size_t> bandCounts;
		for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
		{
			if(bandCounts.count(i->bandname)!=0)
				bandCounts[i->bandname]++;
			else
				bandCounts.insert(std::make_pair(i->bandname, 1));
		}
		size_t maxCount = 0;
		for(std::map<std::string, size_t>::const_iterator i=bandCounts.begin(); i!=bandCounts.end(); ++i)
			maxCount = std::max(maxCount, i->second);
		std::cout << "Largest band has " << maxCount << " sets.\n";
		
		std::ofstream logFile("sedcombine-log.txt");
		
		Model outputModelAll;
		processList(sets, outputModelAll, logFile);
		outputModelAll.Save(outputModelName);
		
		if(maxCount>1)
		{
			for(size_t setPerBand=1; setPerBand!=maxCount; ++setPerBand)
			{
				std::set<std::vector<SEDSet>> permutations;
				for(size_t permutationNr=0; permutationNr!=maxCount*2; ++permutationNr)
				{
					std::vector<SEDSet> permutedSets;
					do {
						makePermutation(setPerBand, maxCount, sets, permutedSets);
					} while(permutations.count(permutedSets) != 0);
					permutations.insert(permutedSets);
					Model outputModel;
					
					processList(permutedSets, outputModel, logFile);
				}
			}
		}
	}
}
