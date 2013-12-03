#include "fitsreader.h"
#include "uvector.h"

#include <iostream>
#include <cmath>

void reportValues(const ao::uvector<double>& samples)
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
	
	std::cout << "N=" << samples.size() << "\tAvg=" << avg << "\tRMS=" << rms << "\tStddev=" << stddev << "\n";
}

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout << "Syntax: imgstats <fitsfile> <trunc value>\n";
	}
	else {
		FitsReader reader(argv[1]);
		const double truncValue = atof(argv[2]);
		
		const size_t width = reader.ImageWidth(), height = reader.ImageHeight();
		ao::uvector<double> image(width * height);
		reader.Read(image.data());
		
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
