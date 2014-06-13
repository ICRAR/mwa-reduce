#ifndef VO_TABLE_PARSER_H
#define VO_TABLE_PARSER_H

#include <functional>
#include <string>

#include <libxml/xmlreader.h>

class VOTableParser
{
public:
	VOTableParser();
	
	~VOTableParser();
	
	void Parse(const std::string& filename);
	
	void SetReadValueCallback(std::function<void(const std::string&, const std::string&)> callback) { _readValueCallback = callback; }
	void SetStartRowCallback(std::function<void()> callback) { _startRowCallback = callback; }
	void SetEndRowCallback(std::function<void()> callback) { _endRowCallback = callback; }
	
private:
	std::function<void(const std::string&, const std::string&)> _readValueCallback;
	std::function<void()> _startRowCallback, _endRowCallback;
	
	xmlTextReaderPtr _reader;
	
	struct FieldStruct
	{
		std::string name;
	};

	void processNode();
	void parseVOTable();
	void parseResource();
	void parseTable();
	void parseField(FieldStruct& field);
	std::string getNextName();
	std::string getNextText();
	std::string getAttribute(const std::string& name);
};

#endif
