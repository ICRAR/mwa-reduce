#include "model.h"
#include "RaDecCoord.h"

#include <cstdlib>
#include <string>
#include <fstream>
#include <stdexcept>

Model::Model(const char* filename)
{
	std::ifstream file(filename);
	
	while(file.good())
	{
		std::string line;
		getline(file, line);
		size_t lineLen = line.size();
		if(lineLen!=0 && !isCommentSymbol(line[0]))
		{
			size_t linePos = 0;
			std::vector<std::string> tokens;
			while(linePos != lineLen)
			{
				// skip leading white space
				while(linePos!=lineLen && isDelimiter(line[linePos])) ++linePos;
				size_t tokenStart = linePos;
				
				// find next white space
				while(linePos!=lineLen && !isDelimiter(line[linePos])) ++linePos;
				size_t tokenEnd = linePos;
				
				if(tokenEnd != tokenStart)
					tokens.push_back(line.substr(tokenStart, tokenEnd-tokenStart));
			}
			
			if(tokens.size() == 8)
			{
				ModelSource source;
				source.SetName(tokens[0]);
				if(tokens[1] != "point")
					throw std::runtime_error("Invalid source type");
				source.SetType(ModelSource::PointSource);
				source.SetPosRA(RaDecCoord::ParseRA(tokens[2].c_str()));
				source.SetPosDec(RaDecCoord::ParseDec(tokens[3].c_str()));
				long double
					fluxA = strtold(tokens[4].c_str(), 0),
					refFreqA = strtold(tokens[5].c_str(), 0)*1000000.0,
					fluxB = strtold(tokens[6].c_str(), 0),
					refFreqB = strtold(tokens[7].c_str(), 0)*1000000.0;
				source.SetRefFreqA(refFreqA);
				source.SetRefFreqB(refFreqB);
				source.Brightness() = SourceStrength<long double>(fluxA, refFreqA, fluxB, refFreqB);
				_sources.push_back(source);
			} else {
				throw std::runtime_error("Invalid token count on line");
			}
		}
	}
}
