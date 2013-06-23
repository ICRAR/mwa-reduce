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
		
		size_t SourceCount() { return _sources.size(); }
		ModelSource &Source(size_t index) { return _sources[index]; }
		const ModelSource &Source(size_t index) const { return _sources[index]; }
		
		const_iterator begin() const { return _sources.begin(); }
		const_iterator end() const { return _sources.end(); }
		
		iterator begin() { return _sources.begin(); }
		iterator end() { return _sources.end(); }
		
		void Optimize();
		
		void AddSource(const ModelSource &source) { _sources.push_back(source); }
		
		void Save(const char *filename);
		
		enum PolarizationType { FullXY };
		
		double TotalFlux(double frequencyStartHz, double frequencyEndHz, size_t polarizationIndex) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->SED().IntegratedFlux(frequencyStartHz, frequencyEndHz, polarizationIndex);
			
			return flux;
		}
		
	private:
		enum PolarizationType _polarizationType;
		std::vector<ModelSource> _sources;
		
		static bool isCommentSymbol(char c) { return c=='#'; }
		static bool isDelimiter(char c) { return c==' ' || c=='\t' || c=='\r' || c=='\n';	}
		void addOptimized(const ModelSource &source);
};

#endif
