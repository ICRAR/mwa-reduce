#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <iomanip>

#include "stopwatch.h"

struct TimingRecord {
	double average, stddev;
	size_t n;
	double sum, sumSq;
	
	void Recalculate()
	{
		average = sum / double(n);
		double sumMeanSquared = sum * sum / n;
		stddev = sqrt((sumSq - sumMeanSquared) / n);
	}
	
	std::string ToLine(const std::string& description, const double val)
	{
		std::ostringstream newLine;
		newLine.precision(11);
		newLine << description << '\t' << val
			<< '\t' << average << '\t' << stddev
			<< '\t' << n << '\t' << sum << '\t' << sumSq << '\n';
		return newLine.str();
	}
};

bool parseLine(const std::string& line, std::string& description, double& value, TimingRecord& record)
{
	std::string tokens[7];
	size_t curPos = 0;
	for(size_t i=0; i!=6; ++i)
	{
		size_t delimitor = line.find('\t', curPos);
		if(delimitor == std::string::npos)
			return false;
		tokens[i] = line.substr(curPos, delimitor - curPos);
		curPos = delimitor+1;
	}
	tokens[6] = line.substr(curPos);
	
	description = tokens[0];
	value = atof(tokens[1].c_str());
	record.average = atof(tokens[2].c_str());
	record.stddev = atof(tokens[3].c_str());
	record.n = atof(tokens[4].c_str());
	record.sum = atof(tokens[5].c_str());
	record.sumSq = atof(tokens[6].c_str());
	return true;
}

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		std::cout << "Syntax: watch <file> <task-description> <task-value> <cmd> [cmd-params ..]\n";
		return -1;
	}
	std::ostringstream str;
	size_t cmdIndex = 4;
	str << '\"' << argv[cmdIndex] << '\"';
	for(int i=cmdIndex+1; i!=argc; ++i)
		str << ' ' << '\"' << argv[i] << '\"';
	const std::string command = str.str();
	const std::string taskDescription(argv[2]);
	const double taskValue(atof(argv[3]));
	
	std::cout << "Executing: " << command << '\n';
	Stopwatch watch(true);
	system(command.c_str());
	double watchTime = watch.Seconds();
	std::cout << "Task \"" << taskDescription << "\" finished in " << watch.ToString() << ".\n";
	
	const char* filename(argv[1]);
	std::ifstream file(filename);
	bool taskFound = false;
	std::map<double, std::string> data;
	while(file.good())
	{
		std::string line;
		std::getline(file, line);
		if(!line.empty())
		{
			std::string curDesc;
			double curVal;
			TimingRecord record;
			if(parseLine(line, curDesc, curVal, record))
			{
				if(!taskFound && curDesc == taskDescription)
				{
					taskFound = true;
					record.n++;
					record.sum += watchTime;
					record.sumSq += watchTime*watchTime;
					record.Recalculate();
					data.insert(std::make_pair(curVal, record.ToLine(taskDescription, curVal)));
				}
				else {
					data.insert(std::make_pair(curVal, line));
				}
			}
		}
	}
	if(!taskFound)
	{
		TimingRecord record;
		record.n = 1;
		record.sum = watchTime;
		record.sumSq = watchTime*watchTime;
		record.average = watchTime;
		record.stddev = std::numeric_limits<double>::quiet_NaN();
		data.insert(std::make_pair(taskValue, record.ToLine(taskDescription, taskValue)));
	}
	file.close();
	
	std::ofstream ofile(filename);
	for(std::map<double, std::string>::const_iterator i=data.begin(); i!=data.end(); ++i)
		ofile << i->second << '\n';
}
