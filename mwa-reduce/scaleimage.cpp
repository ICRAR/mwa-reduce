#include "parser/votableparser.h"
#include "fitsreader.h"
#include "uvector.h"
#include "fitswriter.h"

#include <cmath>
#include <iostream>
#include <limits>

static std::string modelFluxCol, measuredFluxCol;
double curRowModelFlux, curRowMeasuredFlux;
double totalModelFlux, totalMeasuredFlux;
size_t sourcesUsed;

void readValue(const std::string& columnName, const std::string& value)
{
	if(columnName == modelFluxCol)
		curRowModelFlux = atof(value.c_str());
	else if(columnName == measuredFluxCol)
		curRowMeasuredFlux = atof(value.c_str());
}

void startRow()
{
	curRowModelFlux = std::numeric_limits<double>::quiet_NaN();
	curRowMeasuredFlux = std::numeric_limits<double>::quiet_NaN();
}

void endRow()
{
	if(std::isfinite(curRowModelFlux) && std::isfinite(curRowMeasuredFlux))
	{
		totalModelFlux += curRowModelFlux;
		totalMeasuredFlux += curRowMeasuredFlux;
		++sourcesUsed;
	}
}

int main(int argc, char* argv[])
{
	if(argc != 6)
		std::cout << "Syntax: scaleimage <cross-matched vo table> <model flux column name> <measured flux column name> <input fits> <output fits>\n";
	else {
		modelFluxCol = argv[2];
		measuredFluxCol = argv[3];
		VOTableParser parser;
		parser.SetReadValueCallback(&readValue);
		parser.SetStartRowCallback(&startRow);
		parser.SetEndRowCallback(&endRow);
		
		totalModelFlux = 0.0;
		totalMeasuredFlux = 0.0;
		sourcesUsed = 0;
		
		parser.Parse(argv[1]);
		
		std::cout <<
			"Sum over measured   \n"
			"----------------- = " << (totalMeasuredFlux / totalModelFlux) << "\n"
			"  Sum over model\n";
		double factor = totalModelFlux / totalMeasuredFlux;
		
		FitsReader reader(argv[4]);
		ao::uvector<double> values(reader.ImageWidth() * reader.ImageHeight());
		reader.Read(values.data());
		
		for(ao::uvector<double>::iterator i=values.begin(); i!=values.end(); ++i)
			*i *= factor;
		
		FitsWriter writer(reader);
		writer.Write(argv[5], values.data());
	}
	return 0;
}
