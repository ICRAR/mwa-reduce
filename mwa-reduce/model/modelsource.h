#ifndef MODELSOURCE_H
#define MODELSOURCE_H

#include <algorithm>
#include <string>
#include <sstream>

#include "modelcomponent.h"

class ModelSource
{
	public:
		typedef std::vector<ModelComponent>::iterator iterator;
		typedef std::vector<ModelComponent>::const_iterator const_iterator;
		
		ModelSource() : _name(), _components(), _userdata(0)
		{
		}
		
		ModelSource(const ModelSource &source) : _name(source._name), _components(source._components), _userdata(source._userdata), _clusterName(source._clusterName)
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
			_clusterName = source._clusterName;
			return *this;
		}
		
		const std::string &Name() const { return _name; }
		
		void SetName(const std::string &name) { _name = name; }
		
		/**
		 * Returns nullptr in case the source is not part of a cluster.
		 */
		const std::string& ClusterName() const { return _clusterName; }
		
		void SetClusterName(const std::string& clusterName) { _clusterName = clusterName; }
		
		std::string ToString() const;
		
		bool operator<(const ModelSource &rhs) const
		{
			return TotalFlux(Polarization::StokesI)
				< rhs.TotalFlux(Polarization::StokesI);
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
				(*this) += *i;
		}
		
		void operator*=(double factor)
		{
			for(iterator i=begin(); i!=end(); ++i)
				(*i) *= factor;
		}
	
		void CombineMeasurements(const ModelSource& source)
		{
			for(const_iterator i = source.begin(); i!=source.end(); ++i)
				CombineMeasurements(*i);
		}

		void CombineMeasurements(const ModelComponent& component)
		{
			for(iterator i = begin(); i!=end(); ++i)
			{
				if(component.PosDec() == i->PosDec() && component.PosRA() == i->PosRA())
				{
					i->MSED().CombineMeasurements(component.MSED());
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
		
		double TotalFlux(PolarizationEnum polarization) const
		{
			if(_components.empty())
				return 0.0;
			else
				return TotalFlux(begin()->SED().ReferenceFrequencyHz(), polarization);
		}
		
		size_t ComponentCount() const { return _components.size(); }
		
		const class ModelComponent& Component(size_t index) const { return _components[index]; }
		
		void *UserData() const { return _userdata; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		/*void MakeUnitFlux()
		{
			double totalFlux = 0.0;
			double freq = (Peak().MSED().LowestFrequency() + Peak().MSED().HighestFrequency()) * 0.5;
			for(iterator i=begin(); i!=end(); ++i)
			{
				totalFlux += TotalFlux(freq, Polarization::StokesI);
			}
			for(iterator i=begin(); i!=end(); ++i)
			{
				double thisFlux = i->SED().FluxAtFrequency(freq, Polarization::StokesI);
				i->SetSED(MeasuredSED(thisFlux / totalFlux, freq));
			}
		}*/
		
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
				i->SetSED(MeasuredSED(thisFlux * scaleFactor, frequency));
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
					MeasuredSED sed;
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
					MeasuredSED sed;
					sed.AddMeasurement(m);
					i->SetSED(sed);
				}
			}
		}
		
		double MeanRA() const
		{
			std::vector<double> raValues;
			raValues.reserve(_components.size());
			for(const_iterator c=begin(); c!=end(); ++c)
				raValues.push_back(c->PosRA());
			return ImageCoordinates::MeanRA(raValues);
		}
		
		double MeanDec() const
		{
			double sum = 0.0;
			for(const_iterator c=begin(); c!=end(); ++c)
				sum += c->PosDec();
			return sum / _components.size();
		}
		
		void SortComponents()
		{
			std::sort(_components.rbegin(), _components.rend());
		}
		
		bool HasValidMeasurement() const
		{
			for(const_iterator i=begin(); i!=end(); ++i)
				if(i->HasValidMeasurement()) return true;
			return false;
		}
		
		MeasuredSED GetIntegratedMSED() const
		{
			if(_components.empty())
				return MeasuredSED();
			const_iterator i=begin();
			MeasuredSED sum(i->MSED());
			++i;
			while(i != end())
			{
				const MeasuredSED& sed = i->MSED();
				MeasuredSED::const_iterator sedIter = sed.begin();
				MeasuredSED::iterator sumIter = sum.begin();
				while(sedIter != sed.end() && sumIter != sum.end())
				{
					double frequency = sumIter->second.FrequencyHz();
					if(sedIter->second.FrequencyHz() != frequency)
						throw std::runtime_error("GetIntegratedSED() called for source with components having different SED frequency gridding");
					sumIter->second += sedIter->second;
					++sedIter;
					++sumIter;
				}
				++i;
			}
			return sum;
		}
		
	private:
		std::string _name;
		std::vector<ModelComponent> _components;
		void *_userdata;
		std::string _clusterName;
};

class ModelCluster
{
public:
	const std::string& Name() const { return _name; }
	
	void SetName(const std::string& name) { _name = name; }
	
private:
	std::string _name;
};

class SourceGroup
{
public:
		typedef std::vector<ModelSource>::iterator iterator;
		typedef std::vector<ModelSource>::const_iterator const_iterator;
		iterator begin() { return _sources.begin(); }
		iterator end() { return _sources.end(); }
		const_iterator begin() const { return _sources.begin(); }
		const_iterator end() const { return _sources.end(); }
		
		size_t SourceCount() const { return _sources.size(); }
		
		void AddSource(const ModelSource& source) { _sources.push_back(source); }
		
		double TotalFlux(PolarizationEnum polarization) const
		{
			double f = 0.0;
			for(const_iterator s=_sources.begin(); s!=_sources.end(); ++s)
				f += s->TotalFlux(polarization);
			return f;
		}
		
		double MeanRA() const
		{
			std::vector<double> raValues;
			raValues.reserve(_sources.size());
			for(const_iterator s=begin(); s!=end(); ++s)
				raValues.push_back(s->MeanRA());
			return ImageCoordinates::MeanRA(raValues);
		}
		
		double MeanDec() const
		{
			double sum = 0.0;
			for(const_iterator s=_sources.begin(); s!=_sources.end(); ++s)
				sum += s->MeanDec();
			return sum / _sources.size();
		}
private:
	std::vector<ModelSource> _sources;
};

inline std::string ModelSource::ToString() const
{
	std::stringstream s;
	s << "source {\n  name \"" << _name << "\"\n";
	if(!_clusterName.empty())
		s << "  cluster \"" << _clusterName << "\"\n";
	for(const_iterator i=begin(); i!=end(); ++i)
		s << i->ToString();
	s << "}\n";
	return s.str();
}

#endif
