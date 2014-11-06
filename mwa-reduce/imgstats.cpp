#include "fitsreader.h"
#include "uvector.h"

#include <iostream>
#include <cmath>

void reportValues(const ao::uvector<double>& samples, bool useLimit=false, double limitFlux=0.0)
{
	double rms = 0.0;
	double avg = 0.0;
	
	for(ao::uvector<double>::const_iterator i=samples.begin(); i!=samples.end(); ++i)
	{
		rms += *i * *i;
		avg += *i;
	}
	avg /= samples.size();
	rms = sqrt(rms / samples.size());
	
	double stddev = 0.0;
	for(ao::uvector<double>::const_iterator i=samples.begin(); i!=samples.end(); ++i)
	{
		double val = *i - avg;
		stddev += val*val;
	}
	stddev = sqrt(stddev / samples.size());
	
	if(useLimit)
	{
		if(rms > limitFlux)
			std::cout << "1\n";
	}
	else
		std::cout << "N=" << samples.size() << "\tAvg=" << avg << "\tRMS=" << rms << "\tStddev=" << stddev << "\n";
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout << "Syntax: imgstats [-test <rmslimit>] <fitsfile> <trunc value>\n";
	}
	else {
		size_t argi = 1;
		bool hasLimit = false;
		double limitRMS = 0.0;
		
		while(argv[argi][0] == '-')
		{
			std::string p = &argv[argi][1];
			if(p == "test")
			{
				hasLimit = true;
				++argi;
				limitRMS = atof(argv[argi]);
			}
			else throw std::runtime_error("Unknown parameter");
			++argi;
		}
		
		FitsReader reader(argv[argi]);
		const double truncValue = atof(argv[argi+1]);
		
		const size_t width = reader.ImageWidth(), height = reader.ImageHeight();
		ao::uvector<double> image(width * height);
		reader.Read(image.data());
		
		if(hasLimit)
		{
			size_t boxSize = 160;
			ao::uvector<double> samplesAll;
			for(size_t boxY=0; boxY!=boxSize; ++boxY)
			{
				double* ptr = &image[(boxY + (height - boxSize)/2) * width + (width + boxSize)/2];
				for(size_t boxX=0; boxX!=boxSize; ++boxX)
				{
					samplesAll.push_back(*ptr);
					++ptr;
				}
			}
			reportValues(samplesAll, true, limitRMS);
		}
		else {
			size_t boxSize = 20;
			while(boxSize < width && boxSize < height)
			{
				ao::uvector<double> samplesTruncated, samplesAll;
				for(size_t boxY=0; boxY!=boxSize; ++boxY)
				{
					double* ptr = &image[(boxY + (height - boxSize)/2) * width + (width + boxSize)/2];
					for(size_t boxX=0; boxX!=boxSize; ++boxX)
					{
						samplesAll.push_back(*ptr);
						if(*ptr <= truncValue)
							samplesTruncated.push_back(*ptr);
						++ptr;
					}
				}
				std::cout << "Box size " << boxSize << '\n';
				std::cout << "All:       ";
				reportValues(samplesAll);
				std::cout << "Truncated: ";
				reportValues(samplesTruncated);
				boxSize *= 2;
			}
		}
	}
}
