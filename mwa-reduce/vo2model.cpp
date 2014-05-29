#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <libxml/xmlreader.h>
#include <boost/config/posix_features.hpp>
#include "modelsource.h"
#include "model.h"

struct FieldStruct
{
	std::string name;
};

struct ParserData {
	xmlTextReaderPtr reader;
	std::string fluxColumn;
	bool useConstantSI, useSecondFrequency;
	double frequency, secondFrequency, constantSI;
	std::string siColumn, secondFluxCol;
};

std::string getNextName(ParserData& data)
{
	int ret = xmlTextReaderRead(data.reader);
	if (ret != 1)
		throw std::runtime_error("Failed to parse");
		
	const xmlChar *name;
	name = xmlTextReaderConstName(data.reader);
	if (name == NULL)
		name = BAD_CAST "--";
	//std::cout << name << '\n';
	return std::string(reinterpret_cast<const char*>(name));
}

std::string getNextText(ParserData& data)
{
	xmlChar *val = xmlTextReaderReadString(data.reader);
	std::string s(reinterpret_cast<const char*>(val));
	xmlFree(val);
	return s;
}

double getNextTextAsDouble(ParserData& data)
{
	return atof(getNextText(data).c_str());
}

std::string getAttribute(ParserData& data, const std::string& name)
{
	xmlChar* str = xmlTextReaderGetAttribute(data.reader, reinterpret_cast<const xmlChar*>(name.c_str()));
	std::string s(reinterpret_cast<const char*>(str));
	xmlFree(str);
	return s;
}

void parseField(ParserData& data, FieldStruct& field)
{
	bool endElement = false;
	field.name = getAttribute(data, "name");
	do {
		std::string name = getNextName(data);
		if(name == "FIELD" && xmlTextReaderNodeType(data.reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void parseTable(ParserData& data, Model& model)
{
	std::vector<FieldStruct> fields;
	size_t columnIndex = 0;
	bool endElement = false;
	ModelSource source;
	double flux = 0.0, secondFlux = 0.0, spectInd = 0.0;
	
	std::cout << "Table\n";
	do {
		std::string name = getNextName(data);
		if(name == "FIELD") {
			fields.push_back(FieldStruct());
			parseField(data, fields.back());
		}
		else if(name == "TD") {
			if(xmlTextReaderNodeType(data.reader)==15 /*end*/) {
				++columnIndex;
			} else {
				std::string colName = fields[columnIndex].name;
				if(colName == "Name")
					source.SetName(getNextText(data));
				else if(colName == "RAJ2000")
					source.front().SetPosRA(getNextTextAsDouble(data)*(M_PI/180.0));
				else if(colName == "DEJ2000")
					source.front().SetPosDec(getNextTextAsDouble(data)*(M_PI/180.0));
				else if(colName == data.fluxColumn)
					flux = getNextTextAsDouble(data);
				else if(colName == data.siColumn)
					spectInd = getNextTextAsDouble(data);
				else if(data.useSecondFrequency && colName == data.secondFluxCol)
					secondFlux = getNextTextAsDouble(data);
			}
		}
		else if(name == "TR")
		{
			if(xmlTextReaderNodeType(data.reader)==15 /*end*/) {
				if(data.useConstantSI)
					source.front().SED().AddMeasurement(flux, data.frequency, data.constantSI);
				else if(data.useSecondFrequency)
				{
					source.front().SED().AddMeasurement(flux, data.frequency);
					source.front().SED().AddMeasurement(secondFlux, data.secondFrequency);
				}
				else if(spectInd != 0.0)
					source.front().SED().AddMeasurement(flux, data.frequency, spectInd);
				else
					source.front().SED().AddMeasurement(flux, data.frequency);
				model.AddSource(source);
			} else {
				source = ModelSource();
				source.AddComponent(ModelComponent());
				columnIndex = 0;
				spectInd = 0.0;
			}
		}
		else if(name == "TABLE" && xmlTextReaderNodeType(data.reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void parseResource(ParserData& data, Model& model)
{
	bool endElement = false;
	std::cout << "Resource\n";
	do {
		std::string name = getNextName(data);
		if(name == "TABLE")
			parseTable(data, model);
		else if(name == "RESOURCE" && xmlTextReaderNodeType(data.reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void parseVOTable(ParserData& data, Model& model)
{
	bool endElement = false;
	std::cout << "VOTable\n";
	do {
		std::string name = getNextName(data);
		if(name == "RESOURCE")
			parseResource(data, model);
		else if(name == "VOTABLE" && xmlTextReaderNodeType(data.reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void processNode(ParserData& data, Model& model)
{
	std::string name;
	do {
		name = getNextName(data);
	} while(name != "VOTABLE");
	parseVOTable(data, model);
}

void readFile(const std::string& filename, Model& model, ParserData& data)
{
	data.reader = xmlReaderForFile(filename.c_str(), NULL, 0);
	if(data.reader != NULL) {
		processNode(data, model);
		xmlFreeTextReader(data.reader);
	} else
		throw std::runtime_error("Could not open " + filename);
}

int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Syntax: vo2model [-si <val>] [-sicol <SI col>] [-flux2 <colname> <MHz>] <flux col> <frequency in MHz> <in.vot> <out.txt>\n";
		return -1;
	}

	ParserData data;
	data.useConstantSI = false;
	data.useSecondFrequency = false;
	
	size_t argi = 1;
	while(argv[argi][0]=='-')
	{
		std::string p(&argv[argi][1]);
		if(p == "si")
		{
			++argi;
			data.constantSI = atof(argv[argi]);
			data.useConstantSI = true;
		}
		else if(p == "sicol")
		{
			++argi;
			data.siColumn = argv[argi];
		}
		else if(p == "flux2")
		{
			++argi;
			data.secondFluxCol = argv[argi];
			++argi;
			data.secondFrequency = atof(argv[argi])*1e6;
			data.useSecondFrequency = true;
		}
		++argi;
	}
	std::string
		fluxDensityColumn(argv[argi]),
		inFile(argv[argi+2]),
		outFile(argv[argi+3]);
	double
		frequency = atof(argv[argi+1])*1e6;
	Model model;

	LIBXML_TEST_VERSION

	data.fluxColumn = fluxDensityColumn;
	data.frequency = frequency;
	
	readFile(inFile, model, data);

	xmlCleanupParser();
	
	model.Save(outFile.c_str());
	
	return 0;
}
