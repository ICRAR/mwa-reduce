#ifndef MODEL_H
#define MODEL_H

#include <algorithm>
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
			_sources(source._sources),
			_clusters(source._clusters)
		{
		}
		
		Model(const char *filename) { read(filename); }
		Model(const std::string& filename) { read(filename.c_str()); }
		
		void operator=(const Model &source)
		{
			_polarizationType = source._polarizationType;
			_sources = source._sources;
			_clusters = source._clusters;
		}
		
		void operator+=(const Model &rhs);
		
		size_t SourceCount() const { return _sources.size(); }
		
		size_t ClusterCount() const { return _clusters.size(); }
		
		bool Empty() const { return _sources.size() == 0; }
		
		ModelSource &Source(size_t index) { return _sources[index]; }
		const ModelSource &Source(size_t index) const { return _sources[index]; }
		
		const_iterator begin() const { return _sources.begin(); }
		const_iterator end() const { return _sources.end(); }
		
		iterator begin() { return _sources.begin(); }
		iterator end() { return _sources.end(); }
		
		static size_t npos;
		
		void Optimize();
		
		void AddSource(const ModelSource& source) { _sources.push_back(source); }
		
		/**
		 * Returns the index of the source with the given name, or npos if not found.
		 */
		size_t FindSourceIndex(const std::string& sourceName) const
		{
			for(const_iterator i=begin(); i!=end(); ++i)
			{
				if(i->Name() == sourceName)
					return i - begin();
			}
			return npos;
		}
		
		void AddCluster(const ModelCluster& cluster) {
			bool wasAdded = _clusters.insert(std::make_pair(cluster.Name(), cluster)).second;
			if(!wasAdded)
				throw std::runtime_error("AddCluster(): Adding cluster with duplicate name");
		}
		
		ModelCluster& FindOrAddCluster(const std::string& clusterName) {
			std::map<std::string,ModelCluster>::iterator c =
				_clusters.find(clusterName);
			if(c == _clusters.end())
			{
				ModelCluster newCluster;
				newCluster.SetName(clusterName);
				c = _clusters.insert(std::make_pair(clusterName, newCluster)).first;
			}
			return c->second;
		}
		
		ModelCluster& FindCluster(const std::string& clusterName) {
			std::map<std::string,ModelCluster>::iterator c =
				_clusters.find(clusterName);
			if(c == _clusters.end())
				throw std::runtime_error("FindCluster(): Non-existing cluster '" + clusterName + "' requested.");
			return c->second;
		}
		
		const ModelCluster& FindCluster(const std::string& clusterName) const {
			std::map<std::string,ModelCluster>::const_iterator c =
				_clusters.find(clusterName);
			if(c == _clusters.end())
				throw std::runtime_error("FindCluster(): Non-existing cluster '" + clusterName + "' requested.");
			return c->second;
		}
		
		void GetSourcesInCluster(const std::string& clusterName, SourceGroup &destGroup) const
		{
			for(const_iterator s=begin(); s!=end(); ++s)
			{
				if(s->ClusterName() == clusterName)
					destGroup.AddSource(*s);
			}
		}
		
		void GetClusterNames(std::vector<std::string>& clusterNames) {
			clusterNames.clear();
			clusterNames.reserve(_clusters.size());
			std::set<std::string> clusterWasAdded;
			for(const_iterator s=begin(); s!=end(); ++s)
			{
				if(clusterWasAdded.count(s->ClusterName()) == 0)
				{
					clusterNames.push_back(s->ClusterName());
					clusterWasAdded.insert(s->ClusterName());
				}
			}
		}
		
		void RemoveSource(size_t index) { _sources.erase(_sources.begin() + index); }
		
		void Save(const std::string& filename) const { Save(filename.c_str()); }
		void Save(const char* filename) const;
		void Save(std::ostream& stream) const;
		
		enum PolarizationType { FullXY };
		
		double TotalFlux(double frequencyStartHz, double frequencyEndHz, PolarizationEnum polarization) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->TotalFlux(frequencyStartHz, frequencyEndHz, polarization);
			
			return flux;
		}
		
		double TotalFlux(double frequency, PolarizationEnum polarization) const
		{
			double flux = 0.0;
			for(const_iterator i=begin(); i!=end(); ++i)
				flux += i->TotalFlux(frequency, polarization);
			
			return flux;
		}
		
		void SetUnpolarized()
		{
			for(iterator sourcePtr = begin(); sourcePtr!=end(); ++sourcePtr)
			{
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr != sourcePtr->end(); ++compPtr)
				{
					if(compPtr->HasMeasuredSED())
					{
						MeasuredSED& sed = compPtr->MSED();
						for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
						{
							long double totalFlux = m->second.FluxDensity(Polarization::StokesI);
							m->second.SetZeroExceptSinglePol(Polarization::StokesI, totalFlux);
						}
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
		
		void Sort() {
			std::sort(_sources.rbegin(), _sources.rend());
		}
		
		template<class Compare>
		void Sort(Compare comp) {
			std::sort(_sources.rbegin(), _sources.rend(), comp);
		}
	private:
		void read(const char* filename);
		enum PolarizationType _polarizationType;
		std::vector<ModelSource> _sources;
		std::map<std::string, ModelCluster> _clusters;
		
		static bool isCommentSymbol(char c) { return c=='#'; }
		static bool isDelimiter(char c) { return c==' ' || c=='\t' || c=='\r' || c=='\n';	}
		void add(const ModelSource& source);
		void addOptimized(const ModelSource& source);
		void combineMeasurements(const ModelSource& source);
};

#endif
