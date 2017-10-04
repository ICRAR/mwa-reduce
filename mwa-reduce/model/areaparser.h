#ifndef AREA_PARSER_H
#define AREA_PARSER_H

#include "../areaset.h"

#include "../units/radeccoord.h"

#include "tokenizer.h"

class AreaParser : private Tokenizer
{
	public:
		void Parse(AreaSet& areaSet, const std::string& filename)
		{
			std::ifstream str(filename.c_str());
			Parse(areaSet, str);
		}
		
		void Parse(AreaSet& areaSet, std::ifstream& stream)
		{
			SetStream(stream);
			
			std::string line;
			std::getline(stream, line);
			parseVersionLine(line);
			
			std::string token;
			while(getToken(token))
			{
				if(token != "area")
					throw std::runtime_error("Expecting area");
				
				SkyArea area;
				parseArea(area);
				areaSet.AddArea(area);
			}
		}
	private:
		void parseVersionLine(const std::string &line)
		{
		}
		
		void parseArea(SkyArea &area)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "name") area.SetName(getString());
				else if(token == "element") parseElement(area);
				else if(token == "allow-cleaning") area.SetAllowCleaning(getTokenAsBool());
				else if(token == "weight") area.SetWeight(getTokenAsDouble());
				else throw std::runtime_error("Unknown token");
			}
		}
		
		void parseElement(SkyArea &area)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			SkyAreaElement element;
			while(getToken(token) && token != "}")
			{
				if(token == "shape")
				{
					getToken(token);
					if(token == "circle")
					{
						// Read unit
						getToken(token);
						double radius = getTokenAsDouble() / 180.0 * M_PI;
						getToken(token);
						double ra = RaDecCoord::ParseRA(token);
						getToken(token);
						double dec = RaDecCoord::ParseDec(token);
						element.SetToCircle(radius, ra, dec);
					}
					else if(token == "box")
					{
						getToken(token);
						double raUp = RaDecCoord::ParseRA(token);
						getToken(token);
						double decUp = RaDecCoord::ParseDec(token);
						getToken(token);
						double raLo = RaDecCoord::ParseRA(token);
						getToken(token);
						double decLo = RaDecCoord::ParseDec(token);
						element.SetToBox(raUp, raLo, decUp, decLo);
					}
					else throw std::runtime_error("Unknown shape");
				}
				else throw std::runtime_error("Unknown token in area");
			}
			area.AddElement(element);
		}
};

#endif
