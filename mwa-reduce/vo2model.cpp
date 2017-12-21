#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <libxml/xmlreader.h>

#include "model/modelsource.h"
#include "model/model.h"
#include "model/powerlawsed.h"

#include "parser/votableparser.h"

Model model;
ModelSource source;
double flux, spectInd, secondFlux, thirdFlux;

struct FieldStruct
{
	std::string name;
};

struct ParserData {
	xmlTextReaderPtr reader;
	std::string fluxColumn;
	bool useConstantSI, useSecondFrequency, useThirdFrequency, useSIFormat;
	double frequency, secondFrequency, thirdFrequency, constantSI;
	std::string siColumn, secondFluxCol, thirdFluxCol;
} parserData;

double getAsDouble(const std::string& data)
{
	return atof(data.c_str());
}

void readValue(const std::string& columnName, const std::string& value)
{
	if(columnName == "Name")
		source.SetName(value);
	else if(columnName == "RAJ2000" || columnName == "RA_J2000")
		source.front().SetPosRA(getAsDouble(value)*(M_PI/180.0));
	else if(columnName == "DEJ2000" || columnName == "DEC_J2000")
		source.front().SetPosDec(getAsDouble(value)*(M_PI/180.0));
	else if(columnName == parserData.fluxColumn)
		flux = getAsDouble(value);
	else if(columnName == parserData.siColumn)
		spectInd = getAsDouble(value);
	else if(parserData.useSecondFrequency && columnName == parserData.secondFluxCol)
		secondFlux = getAsDouble(value);
	else if(parserData.useThirdFrequency && columnName == parserData.thirdFluxCol)
		thirdFlux = getAsDouble(value);
}

void startRow()
{
	source = ModelSource();
	source.AddComponent(ModelComponent());
	if(!parserData.useSIFormat)
		source.front().SetSED(MeasuredSED());
	flux = std::numeric_limits<double>::quiet_NaN();
	spectInd = std::numeric_limits<double>::quiet_NaN();
	secondFlux = std::numeric_limits<double>::quiet_NaN();
	thirdFlux = std::numeric_limits<double>::quiet_NaN();
}

void endRow()
{
	if(parserData.useConstantSI)
		source.front().MSED().AddMeasurement(flux, parserData.frequency, parserData.constantSI);
	else if(parserData.useSecondFrequency)
	{
		double fluxAlt = secondFlux, freqAlt = parserData.secondFrequency;
		if(!std::isfinite(secondFlux))
		{
			if(!std::isfinite(thirdFlux))
				std::cout << "Warning: flux2 and flux3 unset\n";
			else {
				fluxAlt = thirdFlux;
				freqAlt = parserData.thirdFrequency;
			}
		}
		source.front().MSED().AddMeasurement(flux, parserData.frequency);
		if(std::isfinite(fluxAlt))
			source.front().MSED().AddMeasurement(fluxAlt, freqAlt);
	}
	else if(!std::isfinite(spectInd)) {
		if(parserData.useSIFormat)
		{
			PowerLawSED sed;
			const double fluxes[4] = { flux, 0.0, 0.0, 0.0 };
			ao::uvector<double> slopes{spectInd};
			sed.SetData(parserData.frequency, fluxes, slopes);
			source.front().SetSED(sed);
		}
		else
			source.front().MSED().AddMeasurement(flux, parserData.frequency, spectInd);
	}
	else
		source.front().MSED().AddMeasurement(flux, parserData.frequency);
	model.AddSource(source);
}

void readFile(const std::string& filename)
{
	VOTableParser parser;
	parser.SetReadValueCallback(&readValue);
	parser.SetStartRowCallback(&startRow);
	parser.SetEndRowCallback(&endRow);
	
	parser.Parse(filename);
}

int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Syntax: vo2model [-si <val>] [-sicol <SI col>] [-flux2 <colname> <MHz>] <flux col> <frequency in MHz> <in.vot> <out.txt>\n";
		return -1;
	}

	parserData.useConstantSI = false;
	parserData.useSecondFrequency = false;
	parserData.useSIFormat = false;
	
	size_t argi = 1;
	while(argv[argi][0]=='-')
	{
		std::string p(&argv[argi][1]);
		if(p == "si")
		{
			++argi;
			parserData.constantSI = atof(argv[argi]);
			parserData.useConstantSI = true;
		}
		else if(p == "sicol")
		{
			++argi;
			parserData.siColumn = argv[argi];
			parserData.useSIFormat = true;
		}
		else if(p == "flux2")
		{
			++argi;
			parserData.secondFluxCol = argv[argi];
			++argi;
			parserData.secondFrequency = atof(argv[argi])*1e6;
			parserData.useSecondFrequency = true;
		}
		else if(p == "flux3")
		{
			++argi;
			parserData.thirdFluxCol = argv[argi];
			++argi;
			parserData.thirdFrequency = atof(argv[argi])*1e6;
			parserData.useThirdFrequency = true;
		}
		++argi;
	}
	std::string
		fluxDensityColumn(argv[argi]),
		inFile(argv[argi+2]),
		outFile(argv[argi+3]);
	double
		frequency = atof(argv[argi+1])*1e6;

	LIBXML_TEST_VERSION

	parserData.fluxColumn = fluxDensityColumn;
	parserData.frequency = frequency;
	
	readFile(inFile);

	xmlCleanupParser();
	
	model.Save(outFile.c_str());
	
	return 0;
}
