#include <iostream>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "model/model.h"
#include "model/modelsource.h"
#include "model/powerlawsed.h"

#include "units/radeccoord.h"

struct BBSParseException : public std::runtime_error
{
	BBSParseException(const std::string& s) : std::runtime_error(s)
	{ }
};

struct Headers
{
	int
		nameInd=-1, typeInd=-1, patchInd=-1, raInd=-1, decInd=-1, iInd=-1, spectrInd=-1,
		referenceFrequency=-1,
		majAxisInd=-1, minAxisInd=-1, orientationInd=-1, logSIInd=-1;
};

class BBSLine
{
public:
	BBSLine(const std::string& line) : _startIndex(0), _endIndex(0), _line(line)
	{
	}
	
	bool MoveToNext()
	{
		if(_endIndex >= _line.size())
			return false;
		_startIndex = _endIndex;
		bool insideBracket = false;
		while(_endIndex < _line.size())
		{
			if(_line[_endIndex] == ',' && !insideBracket)
			{
				++_endIndex;
				return true;
			}
			else if(_line[_endIndex] == '[')
				insideBracket = true;
			else if(_line[_endIndex] == ']' && insideBracket)
				insideBracket = false;
			++_endIndex;
		}
		return true;
	}
	
	std::string Value() const {
		if(_endIndex>0)
			return boost::algorithm::trim_copy(_line.substr(_startIndex, _endIndex-_startIndex-1));
		else
			return std::string();
	}
	
private:
	size_t _startIndex, _endIndex;
	std::string _line;
};

void bbs2model(const std::string& input, const std::string& output, const std::string& sourceName)
{
	std::ifstream inFile(input);
	std::string line;
	std::getline(inFile, line);
	if(line.size()>=3 && line.substr(0, 3) == "# (")
		line = line.substr(3);
	else if(boost::to_lower_copy(line.substr(0, 8)) == "format =")
		line = line.substr(8);
	else
		throw BBSParseException("BBS model does not start with format line");
	if(line.size()>=8 && boost::to_lower_copy(line.substr(line.size()-8)) == "= format")
		line = line.substr(0, line.size()-8);
	boost::char_separator<char> sep(" ,()");
	boost::tokenizer<boost::char_separator<char>> tok(line, sep);
	Headers h;
	int index = 0;
	for(auto s : tok)
	{
		std::string key = boost::algorithm::to_lower_copy(boost::algorithm::trim_copy(s));
		std::string defaultVal;
		size_t eq = key.find('=');
		if(eq != std::string::npos)
		{
			if(eq+1 <= key.size())
				defaultVal = key.substr(eq+1);
			key = key.substr(0, eq);
		}
		
		if(key == "patch")
			h.patchInd = index; // skip
		else if(key == "name")
			h.nameInd = index;
		else if(key == "type")
			h.typeInd = index;
		else if(key == "ra")
			h.raInd = index;
		else if(key == "dec")
			h.decInd = index;
		else if(key == "i")
			h.iInd = index;
		else if(key == "spectralterms" || key == "spectralindex")
			h.spectrInd = index;
		else if(key == "referencefrequency")
			h.referenceFrequency = index;
		else if(key == "majoraxis")
			h.majAxisInd = index;
		else if(key == "minoraxis")
			h.minAxisInd = index;
		else if(key == "orientation")
			h.orientationInd = index;
		else if(key == "format" || key == "q" || key == "u" || key == "v")
			; // skip
		else if(key == "logarithmicsi")
			h.logSIInd = index;
		else
			throw BBSParseException("Unknown header: '" + key + "'");
		++index;
	}
	
	Model model;
	std::getline(inFile, line);
	ModelSource globalSource;
	if(!sourceName.empty())
		globalSource.SetName(sourceName);
	while(inFile.good())
	{
		ModelSource source;
		ModelComponent component;
		bool isPatch = false;
		double refFreq = 0;
		ao::uvector<double> frequencyTerms;
		PowerLawSED sed;
		double stokesI = 0.0;
		index = 0;
		
		BBSLine bbsLine(line);
		if(!line.empty() && line[0]!='#')
		{
			while(bbsLine.MoveToNext())
			{
				std::string val = bbsLine.Value();
				if(index == h.nameInd)
					source.SetName(val);
				else if(index == h.patchInd)
					source.SetClusterName(val);
				else if(index == h.raInd)
					component.SetPosRA(RaDecCoord::ParseRA(val));
				else if(index == h.decInd)
					component.SetPosDec(RaDecCoord::ParseDec(val));
				else if(index == h.typeInd)
				{
					std::string typeStr(val);
					boost::algorithm::to_lower(typeStr);
					if(typeStr == "point")
						component.SetType(ModelComponent::PointSource);
					else if (typeStr == "gaussian")
						component.SetType(ModelComponent::GaussianSource);
					else if(typeStr == "")
						isPatch = true;
					else
						throw BBSParseException("Unknown source type: " + val);
				} else if(index == h.spectrInd)
				{
					boost::char_separator<char> freqsep("[,] ");
					boost::tokenizer<boost::char_separator<char>> freqtok(val, freqsep);
					for(auto fval : freqtok)
						frequencyTerms.push_back(atof(fval.c_str()));
				} else if(index == h.majAxisInd)
					component.SetMajorAxis(atof(val.c_str())*(M_PI/180.0/60.0/60.0));
				else if(index == h.minAxisInd)
					component.SetMinorAxis(atof(val.c_str())*(M_PI/180.0/60.0/60.0));
				else if(index == h.orientationInd)
					component.SetPositionAngle(atof(val.c_str())*(M_PI/180.0));
				else if(index == h.referenceFrequency)
					refFreq = atof(val.c_str());
				else if(index == h.iInd)
					stokesI = atof(val.c_str());
				else if(index == h.logSIInd)
					sed.SetIsLogarithmic(val != "false");
				++index;
			}
			double brightness[] = { stokesI, 0.0, 0.0, 0.0 };
			sed.SetData(refFreq, brightness, frequencyTerms);
			component.SetSED(sed);
			if(!isPatch)
			{
				if(sourceName.empty())
				{
					source.AddComponent(component);
					model.AddSource(source);
				}
				else {
					globalSource.AddComponent(component);
				}
			}
		}
		std::getline(inFile, line);
	}
	if(!sourceName.empty())
		model.AddSource(globalSource);
	model.Save(output);
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout <<
			"Syntax: bbs2model [options] <input.txt> <output.txt>\n"
			"Options:\n"
			" -source <name>\n"
			"    Let all components belong to a single source with the given name.\n";
	}
	else {
		int argi = 1;
		std::string sourcename;
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "source")
			{
				++argi;
				sourcename = argv[argi];
			}
			else throw std::runtime_error("Unknown parameter");
			++argi;
		}
		bbs2model(argv[argi], argv[argi+1], sourcename);
	}
	return 0;
}


