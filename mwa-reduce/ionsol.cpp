#include "ionsolutionfile.h"

#include <iostream>

struct Options
{
	bool list;
	bool replaceUnset;
	std::string outputFilename;
	double dlThreshold, dmThreshold, maxGain, minGain;
	bool hasDLThreshold, hasDMThreshold, hasMaxGain, hasMinGain;
};

void flag(IonSolutionFile::Solution& solution)
{
	solution.dl = std::numeric_limits<double>::quiet_NaN();
	solution.dm = std::numeric_limits<double>::quiet_NaN();
	solution.gain = std::numeric_limits<double>::quiet_NaN();
}

void outputSolutionFile(const char* filename, const char* modelFilename, const Options& options)
{
	IonSolutionFile file;
	file.OpenForReading(filename);
	Model model(modelFilename);
	
	std::vector<std::vector<ModelSource*>> sourcesPerDirection;
	file.ReadClusterMetaInfo(model, sourcesPerDirection);

	std::unique_ptr<IonSolutionFile> output;
	if(!options.outputFilename.empty())
	{
		output.reset(new IonSolutionFile());
		output->SetAntennaCount(file.AntennaCount());
		output->SetChannelBlockCount(file.ChannelBlockCount());
		output->SetDirectionCount(file.DirectionCount());
		output->SetEndFrequency(file.EndFrequency());
		output->SetIntervalCount(file.IntervalCount());
		output->SetPolarizationCount(file.PolarizationCount());
		output->SetStartFrequency(file.StartFrequency());
		output->OpenForWriting(options.outputFilename.c_str());
	}
	
	size_t nSources = 0;
	for(size_t i=0; i!=sourcesPerDirection.size(); ++i)
	{
		if(output)
		{
			std::vector<std::string> sourceNames;
			for(size_t s=0; s!=sourcesPerDirection[i].size(); ++s)
				sourceNames.push_back(sourcesPerDirection[i][s]->Name());
			output->WriteClusterMetaInfo(sourcesPerDirection[i][0]->ClusterName(), sourceNames);
		}
		nSources += sourcesPerDirection[i].size();
	}

	if(options.list)
	{
		std::cout
		 << "# nAntenna=" << file.AntennaCount() << '\n'
		 << "# nChannelBlocks=" << file.ChannelBlockCount() << '\n'
		 << "# nDirections=" << file.DirectionCount() << '\n'
		 << "# nIntervals=" << file.IntervalCount() << '\n'
		 << "# nPolarizations=" << file.PolarizationCount() << '\n'
		 << "# nSources=" << nSources << '\n'
		 << "# frequencyStart=" << file.StartFrequency() << '\n'
		 << "# frequencyEnd=" << file.EndFrequency() << '\n'
		 << "#\n# interval channelBlock direction directionname dl dm gain\n";
	}
	
	size_t nFlagged = 0;
	for(size_t interval=0; interval!=file.IntervalCount(); ++interval)
	{
		for(size_t channelBlock=0; channelBlock!=file.ChannelBlockCount(); ++channelBlock)
		{
			for(size_t direction=0; direction!=file.DirectionCount(); ++direction)
			{
				IonSolutionFile::Solution solution;
				file.ReadSolution(solution, interval, channelBlock, 0, direction);
				const std::string& directionName = sourcesPerDirection[direction][0]->ClusterName();
				bool isFlagged = false;
				
				if(options.list)
					std::cout << interval << '\t' << channelBlock << '\t' << direction << '\t' << directionName << '\t' << solution.dl << '\t' << solution.dm << '\t' << solution.gain << '\n';
				
				if(options.replaceUnset)
				{
					if(solution.dl == 0.0 && solution.dm == 0.0 && solution.gain == 1.0)
					{
						flag(solution);
						isFlagged = true;
					}
				}
				
				if(options.hasDLThreshold && std::fabs(solution.dl) > options.dlThreshold)
				{
					flag(solution);
					isFlagged = true;
				}
				
				if(options.hasDMThreshold && std::fabs(solution.dm) > options.dmThreshold)
				{
					flag(solution);
					isFlagged = true;
				}
				
				if(options.hasMinGain && solution.gain < options.minGain)
				{
					flag(solution);
					isFlagged = true;
				}
				
				if(options.hasMaxGain && solution.gain > options.maxGain)
				{
					flag(solution);
					isFlagged = true;
				}
				
				if(output)
					output->WriteSolution(solution, interval, channelBlock, 0, direction);
				
				if(isFlagged)
					++nFlagged;
			}
		}
	}
	
	if(nFlagged!=0)
	{
		std::cout << "# Newly flagged solutions: " << nFlagged << " (" << double(nFlagged) / (file.IntervalCount()*file.ChannelBlockCount()*file.DirectionCount())*100.0 << "% )\n";
	}
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Syntax: ionsol [options] <solutionfile> <model>\n"
			"Options:\n"
			"-list\n"
			"   List all the solutions.\n"
			"-o <name>\n"
			"   Name of output file. Options below have no effect if no output file is given. Output file\n"
			"   should not be the same as the input file.\n"
			"-replace-unset\n"
			"   Replace all dl=0, dm=0, gain=1 values with NaNs.\n"
			"-threshold-dl <val>\n"
			"   Set any value for which abs(dl) > val to NaNs.\n"
			"-threshold-dm <val>\n"
			"   Set any value for which abs(dm) > val to NaNs.\n"
			"-max-gain <val>\n"
			"   Set any value for which gain > val to NaNs.\n"
			"-min-gain <val>\n"
			"   Set any value for which gain < val to NaNs.\n";
	}
	else {
		Options options;
		options.list = false;
		options.replaceUnset = false;
		options.dlThreshold = 0.0; options.hasDLThreshold = false;
		options.dmThreshold = 0.0; options.hasDMThreshold = false;
		options.maxGain = 0.0; options.hasMaxGain = false;
		options.minGain = 0.0; options.hasMinGain = false;
		
		int argi = 1;
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "list")
				options.list = true;
			else if(p == "replace-unset")
				options.replaceUnset = true;
			else if(p == "threshold-dl")
			{
				++argi;
				options.dlThreshold = atof(argv[argi]);
				options.hasDLThreshold = true;
			}
			else if(p == "threshold-dm")
			{
				++argi;
				options.dmThreshold = atof(argv[argi]);
				options.hasDMThreshold = true;
			}
			else if(p == "max-gain")
			{
				++argi;
				options.maxGain = atof(argv[argi]);
				options.hasMaxGain = true;
			}
			else if(p == "min-gain")
			{
				++argi;
				options.minGain = atof(argv[argi]);
				options.hasMinGain = true;
			}
			else if(p == "o")
			{
				++argi;
				options.outputFilename = argv[argi];
			}
			else throw std::runtime_error("Bad parameter");
			++argi;
		}
		
		outputSolutionFile(argv[argi], argv[argi+1], options);
	}
}
