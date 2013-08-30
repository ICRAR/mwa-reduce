#ifndef MODEL_H
#define MODEL_H

#include <vector>

#include "modelsource.h"

class Model
{
	public:
		typedef std::vector<ModelSource>::iterator iterator;
		typedef std::vector<ModelSource>::const_iterator const_iterator;
		
		Model() :
			_polarizationType(FullXY)
		{
		}

		Model(const Model &source) :
			_polarizationType(source._polarizationType),
			_sources(source._sources)
		{
		}
		
		Model(const char *filename);
		
		void operator=(const Model &source)
		{
			_polarizationType = source._polarizationType;
			_sources = source._sources;
		}
		
		void operator+=(const Model &rhs);
		
		size_t SourceCount() const { return _sources.size(); }
		bool Empty() const { return _sources.size() == 0; }
		
		ModelSource &Source(size_t index) { return _sources[index]; }
		const ModelSource &Source(size_t index) const { return _sources[index]; }
		
		const_iterator begin() const { return _sources.begin(); }
		const_iterator end() const { return _sources.end(); }
		
		iterator begin() { return _sources.begin(); }
		iterator end() { return _sources.end(); }
		
		void Optimize();
		
		void AddSource(const ModelSource &source) { _sources.push_back(source); }
		
		void RemoveSource(size_t index) { _sources.erase(_sources.begin() + index); }
		
		void Save(const char* filename);
		void Save(std::ostream& stream);
		
		enum PolarizationType { FullXY };
		
		double TotalFlux(double frequencyStartHz, double frequencyEndHz, size_t polarizationIndex) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->TotalFlux(frequencyStartHz, frequencyEndHz, polarizationIndex);
			
			return flux;
		}
		
		double TotalFlux(double frequency, size_t polarizationIndex) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->TotalFlux(frequency, polarizationIndex);
			
			return flux;
		}
		
		void SortOnBrightness();
		
		void SetUnpolarized()
		{
			for(iterator sourcePtr = begin(); sourcePtr!=end(); ++sourcePtr)
			{
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr != sourcePtr->end(); ++compPtr)
				{
					SpectralEnergyDistribution &sed = compPtr->SED();
					for(SpectralEnergyDistribution::iterator m=sed.begin(); m!=sed.end(); ++m)
					{
						long double totalFlux = 0.0;
						for(size_t p=0; p!=4; ++p)
						{
							long double f = m->second.FluxDensity(p);
							if(std::isfinite(f))
								totalFlux += f;
						}
						totalFlux /= 2.0;
						m->second.SetFluxDensity(0, totalFlux);
						m->second.SetFluxDensity(1, 0.0);
						m->second.SetFluxDensity(2, 0.0);
						m->second.SetFluxDensity(3, totalFlux);
					}
				}
			}
		}
		
		void CombineMeasurements(const Model& model)
		{
			if(Empty())
				(*this) = model;
			else {
				for(const_iterator i = model.begin(); i!=model.end(); ++i)
					combineMeasurements(*i);
			}
		}
		
		size_t ComponentCount() const {
			size_t count = 0;
			for(const_iterator i = begin(); i!=end(); ++i)
				count += i->ComponentCount();
			return count;
		}
	private:
		enum PolarizationType _polarizationType;
		std::vector<ModelSource> _sources;
		
		static bool isCommentSymbol(char c) { return c=='#'; }
		static bool isDelimiter(char c) { return c==' ' || c=='\t' || c=='\r' || c=='\n';	}
		void add(const ModelSource& source);
		void addOptimized(const ModelSource& source);
		void combineMeasurements(const ModelSource& source);
};

#endif
