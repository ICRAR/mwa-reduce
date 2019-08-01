#include "sedanalyser.h"

#include "../../sourcematcher.h"
#include "../../imageadder.h"

#include <iostream>
#include <fstream>
#include <random>
#include <boost/filesystem/operations.hpp>

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
	std::string bandname, shortname, filename, imageFilename;
	double rms, integrationTime;
	Model* model;
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
	
	double totalIntegrationTime = 0.0;
	for(std::vector<SEDSet>::const_iterator setIter=sets.begin(); setIter!=sets.end(); ++setIter)
		totalIntegrationTime += setIter->integrationTime;
	
	//std::vector<Model> bandModels(bandNames.size());
	//std::vector<Model>::iterator bandModelIter = bandModels.begin();
	for(std::set<std::string>::const_iterator bandIter=bandNames.begin(); bandIter!=bandNames.end(); ++bandIter)
	{
		//Model& bandModel = *bandModelIter;
		std::vector<SEDSet> bandSets;
		for(std::vector<SEDSet>::const_iterator setIter=sets.begin(); setIter!=sets.end(); ++setIter)
		{
			if(setIter->bandname == *bandIter)
				bandSets.push_back(*setIter);
		}
		std::cout << "Processing " << bandSets.size() << " sets in band " << *bandIter << '\n';
		
		std::vector<SEDSet>::const_iterator bandSetIter = bandSets.begin();
		Model bandModel = *bandSetIter->model;
		if(bandSetIter->model->Empty())
		{
			std::cout << "Loading " << bandSetIter->shortname << "...\n";
			bandModel = Model(bandSetIter->filename.c_str());
		}
		double weightSum = useRMS ? 1.0/(bandSetIter->rms*bandSetIter->rms) : 1.0;
		
		++bandSetIter;
		while(bandSetIter != bandSets.end())
		{
			double thisWeight = useRMS ? 1.0/(bandSetIter->rms*bandSetIter->rms) : 1.0;
			
			double addInWeight = thisWeight / (weightSum + thisWeight);
			weightSum += thisWeight;
			
			Model outputModel, addedModel = *bandSetIter->model;
			if(addedModel.SourceCount() == 0)
			{
				std::cout << "Loading " << bandSetIter->shortname << "...\n";
				addedModel = Model(bandSetIter->filename.c_str());
			}
			
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
	logFile << totalIntegrationTime << '\t' << sedAnalyser.BestRMS() << '\t' << i->shortname;
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

void ReadList(const std::string& listFilename, std::vector<SEDSet>& sets)
{
	std::ifstream listfile(listFilename);
	
	while(listfile.good())
	{
		SEDSet sedSet;
		listfile >> sedSet.bandname >> sedSet.rms >> sedSet.integrationTime >> sedSet.shortname >> sedSet.filename >> sedSet.imageFilename;
		if(listfile.good())
		{
			sets.push_back(sedSet);
		}
	}
}

size_t MaxSetsInBand(const std::vector<SEDSet>& sets)
{
	// Count nr of sets in the bands
	std::map<std::string, size_t> bandCounts;
	for(std::vector<SEDSet>::const_iterator i=sets.begin(); i!=sets.end(); ++i)
	{
		if(bandCounts.count(i->bandname)!=0)
			bandCounts[i->bandname]++;
		else
			bandCounts.insert(std::make_pair(i->bandname, 1));
	}
	size_t maxSetsInBand = 0;
	for(std::map<std::string, size_t>::const_iterator i=bandCounts.begin(); i!=bandCounts.end(); ++i)
		maxSetsInBand = std::max(maxSetsInBand, i->second);
	std::cout << "Largest band has " << maxSetsInBand << " sets.\n";
	return maxSetsInBand;
}

double measureImagingNoise(const std::vector<SEDSet>& sets, const std::string& filePostfix)
{
	ImageAdder adder;
	for(const SEDSet& s : sets)
	{
		adder.Add(s.imageFilename + filePostfix);
	}
	const std::string filename = "tmp-permutation.fits";
	adder.Finish(filename);
	
	FitsReader reader(filename);
	size_t width = reader.ImageWidth(), height = reader.ImageHeight();
	ao::uvector<double> image(width * height);
	reader.Read(image.data());
	size_t newWidth = std::min<size_t>(1024, width), newHeight = std::min<size_t>(1024, height);
	ao::uvector<double> newImage(newWidth * newHeight);
	for(size_t y=0; y!=newHeight; ++y)
	{
		for(size_t x=0; x!=newWidth; ++x)
		{
			newImage[y*newWidth + x] = image[(y+(height-newHeight)/2)*width + x+(width-newWidth)/2];
		}
	}
	FitsWriter writer(reader);
	writer.SetImageDimensions(newWidth, newHeight);
	writer.Write(filename, newImage.data());
	
	system(("./bane.sh " + filename).c_str());
	
	FitsReader baneReader("bane_out_rms.fits");
	baneReader.Read(newImage.data());
	double sum = 0.0;
	for(double i : newImage)
		sum += i;
	sum /= newWidth*newHeight;
	std::cout << "Noise level: " << sum << '\n';
	
	boost::filesystem::remove(filename);
	boost::filesystem::remove("bane_out_rms.fits");
	return sum;
}

void processImagePermutation(const std::vector<SEDSet>& sets, std::ofstream& outputFile)
{
	double totalIntegrationTime = 0.0;
	for(std::vector<SEDSet>::const_iterator setIter=sets.begin(); setIter!=sets.end(); ++setIter)
		totalIntegrationTime += setIter->integrationTime;
	
	double stokesINoise = measureImagingNoise(sets, "restored-I.fits");
	double stokesVNoise = measureImagingNoise(sets, "image-V.fits");
	outputFile << totalIntegrationTime << '\t' << stokesINoise << '\t' << stokesVNoise << std::endl;
}

void ImageExperiment(const std::string& listFilename, const std::string& outputTextFile)
{
	std::vector<SEDSet> sets;
	ReadList(listFilename, sets);
	
	// Remove non-existing files
	std::vector<SEDSet> purgedList;
	for(SEDSet& s : sets)
	{
		if(boost::filesystem::exists(s.imageFilename + "restored-I.fits"))
			purgedList.push_back(s);
		else
			std::cout << "Not found: " << s.imageFilename << '\n';
	}
	size_t maxSetsInBand = MaxSetsInBand(purgedList);
	
	std::ofstream outputFile(outputTextFile);
	processImagePermutation(purgedList, outputFile);
	
	if(maxSetsInBand>1)
	{
		for(size_t setPerBand=1; setPerBand!=maxSetsInBand; ++setPerBand)
		{
			std::set<std::vector<SEDSet>> permutations;
			for(size_t permutationNr=0; permutationNr!=maxSetsInBand*2; ++permutationNr)
			{
				std::cout << "setPerBand=" << setPerBand << ", permutation " << permutationNr << '\n';
				std::vector<SEDSet> permutedSets;
				do {
					makePermutation(setPerBand, maxSetsInBand, purgedList, permutedSets);
				} while(permutations.count(permutedSets) != 0);
				permutations.insert(permutedSets);
				
				processImagePermutation(permutedSets, outputFile);
			}
		}
	}
}

void SEDExperiment(const std::string& singleSourceName, const std::string& listFilename, const std::string& outputModelName)
{
	std::vector<SEDSet> sets;
	ReadList(listFilename, sets);
	size_t maxSetsInBand = MaxSetsInBand(sets);
	
	for(SEDSet& sedSet : sets)
	{
		std::cout << "Loading " << sedSet.filename << "...\n";
		sedSet.model = new Model(sedSet.filename);
		if(!singleSourceName.empty())
		{
			Model& newModel = *sedSet.model;
			size_t index = newModel.FindSourceIndex(singleSourceName);
			if(index == Model::npos)
				throw std::runtime_error(std::string("Source '") + singleSourceName + "' not found in model");
			ModelSource singleSource = newModel.Source(index);
			newModel = Model();
			newModel.AddSource(singleSource);
		}
	}
			
	std::ofstream logFile("sedcombine-log.txt");
	
	Model outputModelAll;
	processList(sets, outputModelAll, logFile);
	outputModelAll.Save(outputModelName);
	
	if(maxSetsInBand>1)
	{
		for(size_t setPerBand=1; setPerBand!=maxSetsInBand; ++setPerBand)
		{
			std::set<std::vector<SEDSet>> permutations;
			for(size_t permutationNr=0; permutationNr!=maxSetsInBand*2; ++permutationNr)
			{
				std::vector<SEDSet> permutedSets;
				do {
					makePermutation(setPerBand, maxSetsInBand, sets, permutedSets);
				} while(permutations.count(permutedSets) != 0);
				permutations.insert(permutedSets);
				Model outputModel;
				
				processList(permutedSets, outputModel, logFile);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Syntax: sedcombine [-images] [-single <SOURCENAME>] <listfile> <outputmodel>\n\n"
			"format of rows in listfile:\n"
			"bandname RMS integrationTime shortname filename\n";
	}
	else {
		int argi=1;
		std::string singleSourceName;
		bool doProcessImages = false;
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "single")
			{
				++argi;
				singleSourceName = argv[argi];
			}
			else if(p == "images")
				doProcessImages = true;
			++argi;
		}
		if(doProcessImages)
			ImageExperiment(argv[argi], argv[argi+1]);
		else
			SEDExperiment(singleSourceName, argv[argi], argv[argi+1]);
	}
}
