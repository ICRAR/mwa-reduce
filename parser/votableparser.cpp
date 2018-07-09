#include "votableparser.h"

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <libxml/xmlreader.h>

VOTableParser::VOTableParser()
{
	LIBXML_TEST_VERSION
}

VOTableParser::~VOTableParser()
{
	xmlCleanupParser();
}

std::string VOTableParser::getNextName()
{
	int ret = xmlTextReaderRead(_reader);
	if (ret != 1)
		throw std::runtime_error("Failed to parse");
		
	const xmlChar *name;
	name = xmlTextReaderConstName(_reader);
	if (name == NULL)
		name = BAD_CAST "--";
	return std::string(reinterpret_cast<const char*>(name));
}

std::string VOTableParser::getNextText()
{
	xmlChar *val = xmlTextReaderReadString(_reader);
	std::string s;
	if(val != 0)
		s = std::string(reinterpret_cast<const char*>(val));
	xmlFree(val);
	return s;
}

std::string VOTableParser::getAttribute(const std::string& name)
{
	xmlChar* str = xmlTextReaderGetAttribute(_reader, reinterpret_cast<const xmlChar*>(name.c_str()));
	std::string s(reinterpret_cast<const char*>(str));
	xmlFree(str);
	return s;
}

void VOTableParser::parseField(FieldStruct& field)
{
	field.name = getAttribute("name");
	bool endElement = xmlTextReaderIsEmptyElement(_reader);
	while(!endElement) {
		std::string name = getNextName();
		if(name == "FIELD" && xmlTextReaderNodeType(_reader)==15 /*end*/)
			endElement = true;
	}
}

void VOTableParser::parseTable()
{
	std::vector<FieldStruct> fields;
	size_t columnIndex = 0;
	bool endElement = false;
	
	do {
		std::string name = getNextName();
		bool isEmpty = xmlTextReaderIsEmptyElement(_reader);
		if(name == "FIELD") {
			fields.push_back(FieldStruct());
			parseField( fields.back());
		}
		else if(name == "TD") {
			if(xmlTextReaderNodeType(_reader)==15 /*end*/ || isEmpty) {
				++columnIndex;
			} else {
				std::string colName = fields[columnIndex].name;
				std::string value = getNextText();
				_readValueCallback(colName, value);
			}
		}
		else if(name == "TR")
		{
			if(xmlTextReaderNodeType(_reader)==15 /*end*/)
			{
				_endRowCallback();
			}
			else {
				columnIndex = 0;
				_startRowCallback();
			}
		}
		else if(name == "TABLE" && xmlTextReaderNodeType(_reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void VOTableParser::parseResource()
{
	bool endElement = false;
	do {
		std::string name = getNextName();
		if(name == "TABLE")
			parseTable();
		else if(name == "RESOURCE" && xmlTextReaderNodeType(_reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void VOTableParser::parseVOTable()
{
	bool endElement = false;
	do {
		std::string name = getNextName();
		if(name == "RESOURCE")
			parseResource();
		else if(name == "VOTABLE" && xmlTextReaderNodeType(_reader)==15 /*end*/)
			endElement = true;
	} while(!endElement);
}

void VOTableParser::processNode()
{
	std::string name;
	do {
		name = getNextName();
	} while(name != "VOTABLE");
	parseVOTable();
}

void VOTableParser::Parse(const std::string& filename)
{
	_reader = xmlReaderForFile(filename.c_str(), NULL, 0);
	if(_reader != NULL) {
		processNode();
		xmlFreeTextReader(_reader);
	} else
		throw std::runtime_error("Could not open " + filename);
}
