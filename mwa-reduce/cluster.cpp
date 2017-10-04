#include "units/imagecoordinates.h"
#include "uvector.h"

#include "model/model.h"

#include <iostream>
#include <fstream>

class Cluster
{
public:
	bool operator<(const Cluster& rhs) const
	{
		return TotalFlux() < rhs.TotalFlux();
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
			std::vector<double> ras;
			for(std::vector<const ModelSource*>::const_iterator i=_sources.begin(); i!=_sources.end(); ++i)
				ras.push_back((*i)->Peak().PosRA());
			_meanRA = ImageCoordinates::MeanRA(ras);
			_meanDec = 0.0;
			for(std::vector<const ModelSource*>::const_iterator i=_sources.begin(); i!=_sources.end(); ++i)
			{
				const ModelSource& s = **i;
				_meanDec += s.Peak().PosDec();
			}
			_meanDec /= _sources.size();
		}
	}
	
	long double Distance(const ModelSource& source) const
	{
		return AngularDistance(source);
		/*long double dRA = source.Peak().PosRA() - _meanRA, dDec = source.Peak().PosDec() - _meanDec;
		while(dRA > 2.0*M_PI) dRA -= 2.0*M_PI;
		while(dRA < -2.0*M_PI) dRA += 2.0*M_PI;
		if(dRA > M_PI) dRA = 2.0*M_PI - dRA;
		else if(dRA < -M_PI) dRA = 2.0*M_PI + dRA;
		return sqrt(dRA*dRA + dDec*dDec);*/
	}
	long double AngularDistance(const ModelSource& source) const
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
	long double TotalFlux() const
	{
		double sum = 0.0;
		for(std::vector<const ModelSource*>::const_iterator s=_sources.begin(); s!=_sources.end(); ++s)
		{
			double flux = (*s)->TotalFlux(Polarization::StokesI);
			sum += flux;
		}
		return sum;
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
	double sumOfMaxDist = 0.0, maxDistOfAll = 0.0, sumDistSq = 0.0, sumDist = 0.0;
	size_t n = 0;
	for(size_t i=0; i!=clusters.size(); ++i)
	{
		double meanRA = clusters[i].MeanRA();
		if(meanRA < 0.0) meanRA += 2.0 * M_PI;
		std::cout << "Cluster " << (i+1)
			<< " (" << RaDecCoord::RAToString(meanRA) << " " << RaDecCoord::DecToString(clusters[i].MeanDec()) << ")"
			<< ": " << clusters[i].SourceCount() << " sources";
		if(clusters[i].SourceCount() > 0)
		{
			double sum = 0.0;
			double
				maxSource = clusters[i].Source(0)->TotalFlux(Polarization::StokesI),
				minSource = maxSource;
			double maxDist = 0.0;
			for(size_t s=0; s!=clusters[i].SourceCount(); ++s)
			{
				const ModelSource& source = *clusters[i].Source(s);
				double flux = source.TotalFlux(Polarization::StokesI);
				double dist = clusters[i].AngularDistance(source)*180.0/M_PI;
				sumDist += dist;
				sumDistSq += dist*dist;
				if(flux > maxSource) maxSource = flux;
				if(flux < minSource) minSource = flux;
				if(dist > maxDist) maxDist = dist;
				sum += flux;
				++n;
			}
			std::cout << ", min=" << minSource << ", max=" << maxSource << ", sum=" << sum << ", max dist=" << maxDist << " deg";
			sumOfMaxDist += maxDist;
			if(maxDist > maxDistOfAll) maxDistOfAll = maxDist;
		}
		std::cout << '\n';
	}
	std::cout << "Average max dist: " << sumOfMaxDist/clusters.size() << " deg\n";
	std::cout << "Max dist of all : " << maxDistOfAll << " deg\n";
	std::cout << "Average distance: " << sumDist/n << " deg\n";
	std::cout << "Distance RMS    : " << sqrt(sumDistSq/n) << " deg\n";
}

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout << "Syntax: cluster <model-input> <model-output> <clustercount>\n";
		return -1;
	}
	Model model(argv[1]);
	const size_t clusterCount = atoi(argv[3]);
	
	Model prunedModel;
	for(Model::const_iterator i=model.begin(); i!=model.end(); ++i)
	{
		if(i->ComponentCount() != 0)
			prunedModel.AddSource(*i);
	}
	if(model.SourceCount() != prunedModel.SourceCount())
	{
		std::cout << "Removed " << model.SourceCount()-prunedModel.SourceCount() << " sources without components.\n";
	}
	model = prunedModel;
	
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
	
	std::cout << "Angular K-means used " << iterationCount << " iterations.\n";
	
	std::sort(clusters.rbegin(), clusters.rend());
	
	Output(clusters);
	
	Model outputModel;
	for(size_t cIndex=0; cIndex!=clusters.size(); ++cIndex)
	{
		const Cluster& cluster = clusters[cIndex];
		if(cluster.SourceCount() != 0)
		{
			std::ostringstream clusterNameStr;
			clusterNameStr << "cluster" << (cIndex+1);
			ModelCluster modelCluster;
			modelCluster.SetName(clusterNameStr.str());
			outputModel.AddCluster(modelCluster);
			for(size_t s=0; s!=cluster.SourceCount(); ++s)
			{
				const ModelSource& source = *cluster.Source(s);
				ModelSource clusteredSource(source);
				clusteredSource.SetClusterName(modelCluster.Name());
				outputModel.AddSource(clusteredSource);
			}
		}
	}
	outputModel.Save(argv[2]);
	
	std::ofstream str("kvis.ann");
	for(size_t cIndex=0; cIndex!=clusters.size(); ++cIndex)
	{
		const Cluster& cluster = clusters[cIndex];
		str << "TEXT " << cluster.MeanRA()*180.0/M_PI << " " << cluster.MeanDec()*180.0/M_PI << " cluster" << (cIndex+1) << "\n";
	}
	
}
