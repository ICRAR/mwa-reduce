#ifndef MODELSOURCE_H
#define MODELSOURCE_H

#include <algorithm>
#include <string>
#include <sstream>

#include "sourcesdf.h"
#include "sourcesdfwithsamples.h"
#include "radeccoord.h"
#include "spectralenergydistribution.h"

class ModelComponent
{
	public:
		enum Type { PointSource };
		
		ModelComponent() : _type(PointSource), _posRA(0.0), _posDec(0.0), _sed(), _l(0.0), _m(0.0), _userdata(0)
		{
		}
		
		ModelComponent(const ModelComponent& source) : _type(source._type), _posRA(source._posRA), _posDec(source._posDec), _sed(source._sed), _l(source._l), _m(source._m), _userdata(source._userdata)
		{
		}
		
		ModelComponent& operator=(const ModelComponent &source)
		{
			_type = source._type;
			_posRA = source._posRA;
			_posDec = source._posDec;
			_sed = source._sed;
			_l = source._l;
			_m = source._m;
			_userdata = source._userdata;
			return *this;
		}
		
		enum Type Type() const { return _type; }
		long double PosRA() const { return _posRA; }
		long double PosDec() const { return _posDec; }
		const SpectralEnergyDistribution &SED() const { return _sed; }
		long double L() const { return _l; }
		long double M() const { return _m; }
		void *UserData() const { return _userdata; }
		
		void SetType(enum Type type) { _type = type; }
		void SetPosRA(long double posRA) { _posRA = posRA; }
		void SetPosDec(long double posDec) { _posDec = posDec; }
		SpectralEnergyDistribution &SED() { return _sed; }
		void SetSED(const SpectralEnergyDistribution &sed) {
			_sed = sed;
		}
		void SetL(long double l) { _l = l; }
		void SetM(long double m) { _m = m; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		std::string ToString() const
		{
			std::stringstream s;
			s << "  component {\n"
				"    type point\n"
				"    position " << RaDecCoord::RAToString(_posRA) << ' ' << RaDecCoord::DecToString(_posDec) << '\n' <<
				_sed.ToString() << "  }\n";
			return s.str();
		}
		
		bool operator<(const ModelComponent& rhs) const
		{
			return _sed < rhs._sed;
		}
	private:
		enum Type _type;
		long double _posRA, _posDec;
		SpectralEnergyDistribution _sed;
		long double _l, _m;
		void *_userdata;
};

class ModelSource
{
	public:
		typedef std::vector<ModelComponent>::iterator iterator;
		typedef std::vector<ModelComponent>::const_iterator const_iterator;
		
		ModelSource() : _name(), _components(), _userdata(0)
		{
		}
		
		ModelSource(const ModelSource &source) : _name(source._name), _components(source._components), _userdata(source._userdata)
		{
		}
		
		~ModelSource()
		{
		}
		
		ModelSource& operator=(const ModelSource &source)
		{
			_name = source._name;
			_components = source._components;
			_userdata = source._userdata;
			return *this;
		}
		
		const std::string &Name() const { return _name; }
		
		void SetName(const std::string &name) { _name = name; }
		
		std::string ToString() const
		{
			std::stringstream s;
			s << "source {\n  name \"" << _name << "\"\n";
			for(const_iterator i=begin(); i!=end(); ++i)
				s << i->ToString();
			s << "}\n";
			return s.str();
		}
		
		bool operator<(const ModelSource &rhs) const
		{
			return _components[0] < rhs._components[0];
		}
		
		void operator+=(const ModelComponent& rhs)
		{
			for(iterator i = begin(); i!=end(); ++i)
			{
				if(rhs.PosDec() == i->PosDec() && rhs.PosRA() == i->PosRA())
				{
					i->SED() += rhs.SED();
					return;
				}
			}
			_components.push_back(rhs);
		}
	
		void operator+=(const ModelSource& rhs)
		{
			for(const_iterator i = rhs.begin(); i!=rhs.end(); ++i)
			{
				(*this) += *i;
			}
		}
	
		void CombineMeasurements(const ModelSource& source)
		{
			for(const_iterator i = source.begin(); i!=source.end(); ++i)
				combineMeasurements(*i);
		}

		void combineMeasurements(const ModelComponent& component)
		{
			for(iterator i = begin(); i!=end(); ++i)
			{
				if(component.PosDec() == i->PosDec() && component.PosRA() == i->PosRA())
				{
					i->SED().CombineMeasurements(component.SED());
					return;
				}
			}
			throw std::runtime_error("Combining measurements while not same sources were measured!");
		}

		iterator begin() { return _components.begin(); }
		iterator end() { return _components.end(); }
		const_iterator begin() const { return _components.begin(); }
		const_iterator end() const { return _components.end(); }
		
		ModelComponent& front() { return _components.front(); }
		const ModelComponent& front() const { return _components.front(); }
		
		const ModelComponent& Peak() const { return *begin(); }
		ModelComponent& Peak() { return *begin(); }
		
		void AddComponent(const ModelComponent& component) {
			_components.push_back(component);
		}
		
		void ClearComponents() {
			_components.clear();
		}
		
		double TotalFlux(double frequencyStartHz, double frequencyEndHz, PolarizationEnum polarization) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->SED().IntegratedFlux(frequencyStartHz, frequencyEndHz, polarization);
			
			return flux;
		}
		
		double TotalFlux(double frequency, PolarizationEnum polarization) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->SED().FluxAtFrequency(frequency, polarization);
			
			return flux;
		}
		
		size_t ComponentCount() const { return _components.size(); }
		
		void *UserData() const { return _userdata; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		void MakeUnitFlux()
		{
			double totalFlux = 0.0;
			double freq = (Peak().SED().LowestFrequency() + Peak().SED().HighestFrequency()) * 0.5;
			for(iterator i=begin(); i!=end(); ++i)
			{
				totalFlux += TotalFlux(freq, Polarization::StokesI);
			}
			for(iterator i=begin(); i!=end(); ++i)
			{
				double thisFlux = i->SED().FluxAtFrequency(freq, Polarization::StokesI);
				i->SetSED(SpectralEnergyDistribution(thisFlux / totalFlux, freq));
			}
		}
		
		void SetConstantTotalFlux(double newFlux, double frequency)
		{
			double totalFlux = 0.0;
			for(iterator i=begin(); i!=end(); ++i)
			{
				totalFlux += TotalFlux(frequency, Polarization::StokesI);
			}
			double scaleFactor = newFlux / totalFlux;
			for(iterator i=begin(); i!=end(); ++i)
			{
				double thisFlux = i->SED().FluxAtFrequency(frequency, Polarization::StokesI);
				i->SetSED(SpectralEnergyDistribution(thisFlux * scaleFactor, frequency));
			}
		}
		
		void SetConstantTotalFlux(const double* newFluxes, double frequency)
		{
			double totalFlux = fabs(TotalFlux(frequency, Polarization::StokesI));
			
			if(totalFlux == 0.0)
			{
				for(iterator i=begin(); i!=end(); ++i)
				{
					Measurement m;
					m.SetFrequencyHz(frequency);
					for(size_t p=0; p!=4; ++p)
						m.SetFluxDensityFromIndex(p, newFluxes[p] / (double) ComponentCount());
					SpectralEnergyDistribution sed;
					sed.AddMeasurement(m);
					i->SetSED(sed);
				}
			}
			else {
				totalFlux *= 0.5;
				double scaleFactor[4];
			
				for(size_t p=0; p!=4; ++p)
					scaleFactor[p] = newFluxes[p] / totalFlux;
				
				for(iterator i=begin(); i!=end(); ++i)
				{
					Measurement m;
					m.SetFrequencyHz(frequency);
					double thisFlux = 0.5*(i->SED().FluxAtFrequency(frequency, Polarization::StokesI));
					for(size_t p=0; p!=4; ++p)
					{
						m.SetFluxDensityFromIndex(p, thisFlux * scaleFactor[p]);
					}
					SpectralEnergyDistribution sed;
					sed.AddMeasurement(m);
					i->SetSED(sed);
				}
			}
		}
		
		double MeanRA() const
		{
			long double firstRA = _components[0].PosRA();
			bool doWrap =
				(firstRA > - 0.25*M_PI && firstRA < 0.25*M_PI) ||
				(firstRA > 1.75*M_PI && firstRA < 2.25*M_PI) ||
				(firstRA > -2.25*M_PI && firstRA < -1.75*M_PI);

			double sum = 0.0;
			for(const_iterator c=_components.begin(); c!=_components.end(); ++c)
			{
				double ra = c->PosRA();
				if(doWrap)
				{
					if(ra > M_PI) ra -= 2.0 * M_PI;
					else if(ra < -M_PI) ra += 2.0 * M_PI;
				}
				sum += ra;
			}
			return sum / _components.size();
		}
		
		double MeanDec() const
		{
			double sum = 0.0;
			for(const_iterator c=_components.begin(); c!=_components.end(); ++c)
				sum += c->PosDec();
			return sum / _components.size();
		}
		
		void SortComponents()
		{
			std::sort(_components.rbegin(), _components.rend());
		}
	private:
		std::string _name;
		std::vector<ModelComponent> _components;
		void *_userdata;
};

#endif
