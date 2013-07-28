#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include "../model.h"

#include "tokenizer.h"

#include <cstdlib>
#include <fstream>
#include <stdexcept>

class ModelParser : private Tokenizer
{
	public:
		ModelParser()
		{
		}
		
		void Parse(Model &model, std::ifstream &stream)
		{
			SetStream(stream);
			
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
};

#endif

