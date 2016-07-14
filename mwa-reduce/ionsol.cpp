#include "ionsolutionfile.h"

#include <iostream>

void outputSolutionFile(const char* filename, const char* modelFilename)
{
	IonSolutionFile file;
	file.OpenForReading(filename);
	Model model(modelFilename);
	std::cout << "# nAntenna=" << file.AntennaCount() << '\n';
	std::cout << "# nChannelBlocks=" << file.ChannelBlockCount() << '\n';
	std::cout << "# nDirections=" << file.DirectionCount() << '\n';
	std::cout << "# nIntervals=" << file.IntervalCount() << '\n';
	std::cout << "# nPolarizations=" << file.PolarizationCount() << '\n';
	
	std::vector<std::vector<ModelSource*>> sourcesPerDirection;
	file.ReadClusterMetaInfo(model, sourcesPerDirection);

	std::cout << "#\n# interval channelBlock direction directionname dl dm gain\n";
	for(size_t interval=0; interval!=file.IntervalCount(); ++interval)
	{
		for(size_t channelBlock=0; channelBlock!=file.ChannelBlockCount(); ++channelBlock)
		{
			for(size_t direction=0; direction!=file.DirectionCount(); ++direction)
			{
				IonSolutionFile::Solution solution;
				file.ReadSolution(solution, interval, channelBlock, 0, direction);
				const std::string& directionName = sourcesPerDirection[direction][0]->ClusterName();
				std::cout << interval << '\t' << channelBlock << '\t' << direction << '\t' << directionName << '\t' << solution.dl << '\t' << solution.dm << '\t' << solution.gain << '\n';
			}
		}
	}
}

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout << "Syntax: ionsol <solutionfile> <model>\n";
	}
	else {
		outputSolutionFile(argv[1], argv[2]);
	}
}
