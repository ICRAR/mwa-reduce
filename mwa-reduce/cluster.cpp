#include "model.h"
#include "imagecoordinates.h"
#include "uvector.h"

#include <iostream>

class Cluster
{
public:
	Cluster() { }
	Cluster(const Cluster& source) :
		_sources(source._sources),
		_meanRA(source._meanRA),
		_meanDec(source._meanDec)
	{ }
	void operator=(const Cluster& rhs)
	{
		_sources = rhs._sources;
		_meanRA = rhs._meanRA;
		_meanDec = rhs._meanDec;
	}
	
	void AddSource(const ModelSource* source)
	{
		_sources.push_back(source);
	}
	void RemoveSource(const ModelSource* source)
	{
		for(std::vector<const ModelSource*>::iterator i=_sources.begin(); i!=_sources.end(); ++i)
		{
			if(*i == source)
			{
				_sources.erase(i);
				return;
			}
		}
		throw std::runtime_error("Could not find source to be removed.");
	}
	void Clear() { _sources.clear(); }
	size_t SourceCount() const { return _sources.size(); }
	const ModelSource* Source(size_t i) const { return _sources[i]; }
	
	void RecalculateMean()
	{
		if(!_sources.empty())
		{
			_meanRA = 0.0;
			_meanDec = 0.0;
			for(std::vector<const ModelSource*>::const_iterator i=_sources.begin(); i!=_sources.end(); ++i)
			{
				const ModelSource& s = **i;
				_meanRA += s.Peak().PosRA();
				_meanDec += s.Peak().PosDec();
			}
			_meanRA /= _sources.size();
			_meanDec /= _sources.size();
		}
	}
	
	double Distance(const ModelSource& source) const
	{
		double dRA = source.Peak().PosRA() - _meanRA, dDec = source.Peak().PosDec() - _meanDec;
		return sqrt(dRA*dRA + dDec*dDec);
	}
	double AngularDistance(const ModelSource& source) const
	{
		return ImageCoordinates::AngularDistance<long double>(source.Peak().PosRA(), source.Peak().PosDec(), _meanRA, _meanDec);
	}
	long double MeanRA() const { return _meanRA; }
	long double MeanDec() const { return _meanDec; }
	void SetMean(const ModelSource& source)
	{
		_meanRA = source.Peak().PosRA();
		_meanDec = source.Peak().PosDec();
	}
private:
	std::vector<const ModelSource*> _sources;
	
	long double _meanRA, _meanDec;
};

size_t NearestCluster(const std::vector<Cluster>& clusters, const ModelSource& source)
{
	size_t index = 0;
	double minDistance = clusters.front().Distance(source);
	for(size_t i=1; i!=clusters.size(); ++i)
	{
		double d = clusters[i].Distance(source);
		if(d < minDistance)
		{
			index = i;
			minDistance = d;
		}
	}
	return index;
}

void Output(const std::vector<Cluster>& clusters)
{
	for(size_t i=0; i!=clusters.size(); ++i)
	{
		std::cout << "Cluster " << i
			<< " (" << RaDecCoord::RAToString(clusters[i].MeanRA()) << " " << RaDecCoord::DecToString(clusters[i].MeanDec()) << ")"
			<< ": " << clusters[i].SourceCount() << " sources";
		if(clusters[i].SourceCount() > 0)
		{
			double sum = clusters[i].Source(0)->TotalFlux(150000000.0, 0);
			double maxSource = sum, minSource = sum;
			double maxDist = clusters[i].AngularDistance(*clusters[i].Source(0));
			for(size_t s=1; s!=clusters[i].SourceCount(); ++s)
			{
				const ModelSource& source = *clusters[i].Source(s);
				double flux = source.TotalFlux(150000000.0, 0);
				double dist = clusters[i].AngularDistance(source);
				if(flux > maxSource) maxSource = flux;
				if(flux < minSource) minSource = flux;
				if(dist > maxDist) maxDist = dist;
			}
			std::cout << ", min=" << minSource << ", max=" << maxSource << ", sum=" << sum << ", max dist=" << (maxDist*180.0/M_PI);
		}
		std::cout << '\n';
	}
}

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout << "Syntax: cluster <model-input> <model-output> <clustercount>\n";
	}
	const Model& model(argv[1]);
	const size_t clusterCount = atoi(argv[3]);
	
	std::vector<Cluster> clusters(clusterCount);
	size_t clusterIndex = 0;
	ao::uvector<size_t> clusterIndexPerSource(model.SourceCount());
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const ModelSource& s = model.Source(i);
		if(clusters[clusterIndex].SourceCount() == 0)
			clusters[clusterIndex].SetMean(s);
		clusters[clusterIndex].AddSource(&s);
		clusterIndexPerSource[i] = clusterIndex;
		clusterIndex = (clusterIndex + 1) % clusterCount;
	}
	
	bool change;
	size_t iterationCount = 0;
	do {
		change = false;
		
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			size_t nearest = NearestCluster(clusters, model.Source(i));
			if(clusterIndexPerSource[i] != nearest)
			{
				change = true;
				clusters[clusterIndexPerSource[i]].RemoveSource(&model.Source(i));
				clusters[nearest].AddSource(&model.Source(i));
				clusterIndexPerSource[i] = nearest;
			}
		}
		
		for(std::vector<Cluster>::iterator c=clusters.begin(); c!=clusters.end(); ++c)
			c->RecalculateMean();
		
		++iterationCount;
	} while(change && iterationCount < 1000);
	
	std::cout << "K-means used " << iterationCount << " iterations.\n";
	
	Output(clusters);
}
