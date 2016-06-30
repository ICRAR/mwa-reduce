#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include "model.h"
#include "modelsource.h"

#include "tokenizer.h"
#include "powerlawsed.h"

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
			if(stream.bad())
					throw std::runtime_error("Error parsing model");
			
			std::string token;
			while(getToken(token))
			{
				if(token != "source")
					throw std::runtime_error("Expecting source");
				
				ModelSource source;
				parseSource(source);
				model.AddSource(source);
				model.FindOrAddCluster(source.ClusterName());
			}
		}
		
	private:
		bool _fileVersion1_0;
		
		void parseVersionLine(const std::string &line)
		{
			const std::string headerStart = "skymodel fileformat ";
			if(line.substr(0, headerStart.size()) != headerStart)
				throw std::runtime_error("The model file didn't start with a skymodel header");
			std::string version = line.substr(headerStart.size());
			if(version == "1.0")
			{
				_fileVersion1_0 = true;
				std::cout << "Warning: this model is written in the old sky-model file format. The old format does \n"
				             "         not properly support polarization, so is deprecated.\n"
										 "         Use \"editmodel -m new_model.txt old_model.txt\" to convert the model.\n";
			}
			else if(version != "1.1")
			{
				throw std::runtime_error("This model specified file version \"" + version + "\": don't know how to read.\n");
			}
			else _fileVersion1_0 = false;
		}
		
		void parseSource(ModelSource& source)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "name") source.SetName(getString());
				else if(token == "cluster") source.SetClusterName(getString());
				else if(token == "component") {
					ModelComponent component;
					parseComponent(component);
					source.AddComponent(component);
				}
				else throw std::runtime_error("Unknown token");
			}
		}
		
		void parseComponent(ModelComponent &component)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			while(getToken(token) && token != "}")
			{
				if(token == "type") {
					getToken(token);
					if(token == "point")
						component.SetType(ModelComponent::PointSource);
					else if(token == "gaussian")
						component.SetType(ModelComponent::GaussianSource);
					else
						throw std::runtime_error("Unsupported component type");
				}
				else if(token == "position")
				{
					getToken(token);
					component.SetPosRA(RaDecCoord::ParseRA(token));
					getToken(token);
					component.SetPosDec(RaDecCoord::ParseDec(token));
				}
				else if(token == "measurement") {
					Measurement measurement;
					parseMeasurement(measurement);
					if(component.HasMeasuredSED())
						component.MSED().AddMeasurement(measurement);
					else if(component.HasSED())
						throw std::runtime_error("Invalid 'measurement' combined with other brightness specification");
					else {
						component.SetSED(MeasuredSED());
						component.MSED().AddMeasurement(measurement);
					}
				}
				else if(token == "sed") {
					PowerLawSED plSED;
					parsePowerLawSED(plSED);
					if(component.HasSED())
						throw std::runtime_error("Invalid 'sed' combined with other brightness specification");
					else 
						component.SetSED(plSED);
				}
				else if(token == "shape") {
					getToken(token);
					component.SetMajorAxis(atof(token.c_str()) * M_PI / 60.0/60.0/180.0);
					getToken(token);
					component.SetMinorAxis(atof(token.c_str()) * M_PI / 60.0/60.0/180.0);
					getToken(token);
					component.SetPositionAngle(atof(token.c_str()) * M_PI / 180.0);
				}
				else if(token == "major-axis") {
					getToken(token);
					component.SetMajorAxis(atof(token.c_str()) * M_PI / 60.0/60.0/180.0);
				}
				else if(token == "minor-axis") {
					getToken(token);
					component.SetMinorAxis(atof(token.c_str()) * M_PI / 60.0/60.0/180.0);
				}
				else if(token == "position-angle") {
					getToken(token);
					component.SetPositionAngle(atof(token.c_str()) * 180.0/M_PI);
				}
				else throw std::runtime_error("Unknown keyname in component");
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
					if(_fileVersion1_0)
					{
						std::complex<double> linFluxes[4];
						for(size_t p=0; p!=4; ++p)
							linFluxes[p] = getTokenAsDouble();
						double stokesFluxes[4];
						Polarization::LinearToStokes(linFluxes, stokesFluxes);
						for(size_t p=0; p!=4; ++p)
							measurement.SetFluxDensityFromIndex(p, stokesFluxes[p]);
					}
					else {
						for(size_t p=0; p!=4; ++p)
							measurement.SetFluxDensityFromIndex(p, getTokenAsDouble());
					}
				}
				else if(token == "type") {
					getToken(token);
					if(token == "absolute") ; // ignore
					else if(token == "apparent")
						throw std::runtime_error("Model no longer allows apparent values");
					else throw std::runtime_error("Measurement type should be absolute or apparent");
				}
				else if(token == "bandwidth") {
					measurement.SetBandWidthHz(getTokenAsDouble());
					getToken(token);
				}
				else if(token == "fluxdensity-stddev") {
					getToken(token);
					for(size_t p=0; p!=4; ++p)
						measurement.SetFluxDensityStddevFromIndex(p, getTokenAsDouble());
				}
				else if(token == "beam-value") {
					// ignore
				}
				else throw std::runtime_error("Unknown token");
			}
		}
		
		void parsePowerLawSED(PowerLawSED& sed)
		{
			std::string token;
			getToken(token);
			if(token != "{")
				throw std::runtime_error("Expecting {");
			double refFrequency = 0.0;
			double brightness[4] = { 0.0, 0.0, 0.0, 0.0 };
			ao::uvector<double> terms;
			bool hasFrequency = false, hasBrightness = false;
			while(getToken(token) && token != "}")
			{
				if(token == "frequency") {
					if(hasFrequency)
						throw std::runtime_error("Double frequency specification");
					refFrequency = getTokenAsDouble()*1000000.0;
					getToken(token); // unit
					hasFrequency = true;
				}
				else if(token == "fluxdensity") {
					if(hasBrightness)
						throw std::runtime_error("Double brightness specification");
					getToken(token); // unit
					for(size_t p=0; p!=4; ++p)
						brightness[p] = getTokenAsDouble();
					hasBrightness = true;
				}
				else if(token == "spectral-index") {
					if(!terms.empty())
						throw std::runtime_error("Double SI specification");
					getToken(token);
					if(token != "{")
							throw std::runtime_error("Expecting {");
					while(getToken(token) && token != "}")
					{
						terms.push_back(atof(token.c_str()));
					}
				}
				else throw std::runtime_error("Unknown token");
			}
			if(!hasFrequency || !hasBrightness || terms.empty())
				throw std::runtime_error("Incomplete SED specification");
			sed.SetData(refFrequency, brightness, terms);
		}
};

#endif

