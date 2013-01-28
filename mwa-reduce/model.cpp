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
			
			if(tokens.size() >= 8)
			{
				ModelSource source;
				source.SetName(tokens[0]);
				if(tokens[1] != "point")
					throw std::runtime_error("Invalid source type");
				source.SetType(ModelSource::PointSource);
				source.SetPosRA(RaDecCoord::ParseRA(tokens[2].c_str()));
				source.SetPosDec(RaDecCoord::ParseDec(tokens[3].c_str()));
				SourceSDF<long double> *sdf = SourceSDF<long double>::ParseLine(tokens.begin()+4, tokens.end());
				source.SetBrightness(*sdf);
				delete sdf;
				_sources.push_back(source);
			} else {
				throw std::runtime_error("Invalid token count on line");
			}
		}
	}
}

void Model::Optimize()
{
	Model copy(*this);
	_sources.clear();
	for(const_iterator i = copy.begin(); i!=copy.end(); ++i)
		addOptimized(*i);
}

void Model::addOptimized(const ModelSource& source)
{
	for(iterator i = begin(); i!=end(); ++i)
	{
		if(source.PosDec() == i->PosDec() && source.PosRA() == i->PosRA())
		{
			const SourceSDFWithSamples<long double> *sdfSource = dynamic_cast<const SourceSDFWithSamples<long double> *>(&source.Brightness());
			SourceSDFWithSamples<long double> *sdfDest = dynamic_cast<SourceSDFWithSamples<long double> *>(&i->Brightness());
			if(sdfSource == 0 || sdfDest == 0)
				throw std::runtime_error("Can only optimize models with sampled sdfs");
			
			SourceSDFWithSamples<long double>::iterator g = sdfDest->begin();
			for(SourceSDFWithSamples<long double>::const_iterator f = sdfSource->begin();
					f!=sdfSource->end(); ++f, ++g)
			{
				g->second += f->second;
			}
			return;
		}
	}
	_sources.push_back(source);
}
