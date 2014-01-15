#include <fstream>
#include <iostream>
#include <sstream>

#include "stopwatch.h"

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout << "Syntax: watch <file> <task-description> <cmd> [cmd-params ..]\n";
		return -1;
	}
	std::ostringstream str;
	str << argv[3];
	for(int i=4; i!=argc; ++i)
		str << ' ' << argv[i];
	const std::string command = str.str();
	
	std::cout << "Executing: " << command << '\n';
	Stopwatch watch(true);
	system(command.c_str());
	const std::string taskDescription(argv[2]);
	std::cout << "Task \"" << taskDescription << "\" finished in " << watch.ToString() << ".\n";
	
	std::ostringstream newLine;
	newLine << taskDescription << '\t' << watch.Seconds();
	
	const char* filename(argv[1]);
	std::ifstream file(filename);
	bool taskFound = false;
	std::vector<std::string> data;
	while(file.good())
	{
		std::string line;
		std::getline(file, line);
		if(!line.empty())
		{
			data.push_back(line);
			if(!taskFound)
			{
				size_t descEnd = line.find('\t');
				if(descEnd != std::string::npos)
				{
					std::string curDesc = line.substr(0, descEnd);
					if(curDesc == taskDescription)
					{
						taskFound = true;
						data.back() = newLine.str();
					}
				}
			}
		}
	}
	if(!taskFound)
		data.push_back(newLine.str());
	file.close();
	
	std::ofstream ofile(filename);
	for(std::vector<std::string>::const_iterator i=data.begin(); i!=data.end(); ++i)
		ofile << *i << '\n';
}
