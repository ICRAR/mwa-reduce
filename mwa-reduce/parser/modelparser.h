#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include "../model.h"

#include <cstdlib>
#include <fstream>
#include <stdexcept>

class ModelParser
{
	public:
		ModelParser()
		{
		}
		
		void Parse(Model &model, std::ifstream &stream)
		{
			_stream = &stream;
			
			std::string line;
			std::getline(stream, line);
			parseVersionLine(line);
			
			std::string token;
			while(getToken(token))
			{
				if(token != "source")
					throw std::runtime_error("Expecting source");
				
				ModelSource source;
				parseSource(source);
				model.AddSource(source);
			}
		}
		
	private:
		void parseVersionLine(const std::string &line)
		{
		}
		void parseSource(ModelSource &source)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "name") source.SetName(getString());
				else if(token == "component") parseComponent(source);
				else throw std::runtime_error("Unknown token");
			}
		}
		void parseComponent(ModelSource &source)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "type") {
					getToken(token);
					if(token != "point")
						throw std::runtime_error("Unsupported component type");
				}
				else if(token == "position")
				{
					getToken(token);
					source.SetPosRA(RaDecCoord::ParseRA(token));
					getToken(token);
					source.SetPosDec(RaDecCoord::ParseDec(token));
				}
				else if(token == "measurement") {
					Measurement measurement;
					parseMeasurement(measurement);
					source.SED().AddMeasurement(measurement);
				}
			}
		}
		void parseMeasurement(Measurement &measurement)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "frequency") {
					measurement.SetFrequencyHz(getTokenAsDouble()*1000000.0);
					getToken(token);
				}
				else if(token == "fluxdensity") {
					getToken(token); // unit
					for(size_t p=0; p!=4; ++p)
						measurement.SetFluxDensity(p, getTokenAsDouble());
				}
				else if(token == "type") {
					getToken(token);
					if(token == "absolute") measurement.SetIsApparent(false);
					else if(token == "apparent") measurement.SetIsApparent(true);
					else throw std::runtime_error("Measurement type should be absolute or apparent");
				}
				else if(token == "bandwidth") {
					measurement.SetBandWidthHz(getTokenAsDouble());
					getToken(token);
				}
				else if(token == "fluxdensity-stddev") {
					getToken(token);
					for(size_t p=0; p!=4; ++p)
						measurement.SetFluxDensityStddev(p, getTokenAsDouble());
				}
				else if(token == "beam-value") {
					for(size_t p=0; p!=4; ++p)
						measurement.SetBeamValue(p, getTokenAsDouble());
				}
				else throw std::runtime_error("Unknown token");
			}
		}
		std::string getString()
		{
			std::string token;
			getToken(token);
			if(token.size()<2 || token[0]!='\"' || token[token.size()-1]!='\"')
				throw std::runtime_error("Expecting string");
			return token.substr(1, token.size()-2);
		}
		double getTokenAsDouble()
		{
			std::string token;
			getToken(token);
			return atof(token.c_str());
		}
		
		bool getToken(std::string &token)
		{
			if(!_stream->good()) return false;
			
			std::stringstream s;
			bool finished = false;
			enum { StateStart, StateInToken, StateInString, StateEscaped, StateInLineComment } state = StateStart;
			do {
				char c;
				_stream->get(c);
				if(!_stream->fail())
				{
					switch(state)
					{
						case StateStart:
							if(!(c == ' ' || c == '\t' || c == '\n'))
							{
								if(c == '/')
									state = StateEscaped;
								else
								{
									s << c;
									if(c == '\"')
										state = StateInString;
									else
										state = StateInToken;
								}
							}
							break;
						case StateInString:
							s << c;
							if(c == '\"')
								finished = true;
							break;
						case StateInToken:
							if(c == ' ' || c == '\t' || c == '\n' || c == ';')
								finished = true;
							else
								s << c;
							break;
						case StateEscaped:
							if(c == '/')
								state = StateInLineComment;
							else
								throw std::runtime_error("Incorrect /");
						case StateInLineComment:
							if(c == '\n')
								state = StateStart;
							break;
					}
				}
				else {
					token = s.str();
					return !token.empty();
				}
			} while(!finished);
			token = s.str();
			return true;
		}
		
		std::ifstream *_stream;
};

#endif

