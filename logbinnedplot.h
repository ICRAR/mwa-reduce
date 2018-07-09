#ifndef LOG_BINNED_PLOT_H
#define LOG_BINNED_PLOT_H

#include "uvector.h"
#include "gnuplot.h"

#include <cmath>
#include <map>
#include <limits>

class LogBinnedPlot
{
public:
	explicit LogBinnedPlot(size_t binsPerOrderOfMagnitude = 4) : _binsPerOrderOfMagnitude(binsPerOrderOfMagnitude)
	{ }
	
	void Add(double x, double y, double yVariance)
	{
		Bin& bin = _bins[valueToBinCenter(x)];
		double val = y / yVariance;
		bin.sum += val;
		bin.sumSq += val*val;
		bin.weight += 1.0 / yVariance;
	}
	
	void Add(const ao::uvector<double>& x, const ao::uvector<double>& y)
	{
		for(size_t i=0; i!=x.size(); ++i)
		{
			Add(x[i], y[i], 1.0);
		}
	}
	
	void Plot(GNUPlot::Line& line) const
	{
		for(std::map<int, Bin>::const_iterator i=_bins.begin();
				i!=_bins.end(); ++i)
		{
			double avg = i->second.average();
			double sigma = i->second.stdErr();
			line.AddPoint(binToValue(i->first), avg, sigma);
		}
	}
	
	void PlotWithOtherErrorBars(GNUPlot::Line& line, const LogBinnedPlot& plotForErrorBars) const
	{
		std::map<int, Bin>::const_iterator e = plotForErrorBars._bins.begin();
		for(std::map<int, Bin>::const_iterator i=_bins.begin();
				i!=_bins.end(); ++i, ++e)
		{
			double avg = i->second.average();
			double sigma = e->second.stdErr();
			line.AddPoint(binToValue(i->first), avg, sigma);
		}
	}
	
	void GetData(ao::uvector<double>& x, ao::uvector<double>& y)
	{
		x.clear();
		y.clear();
		x.reserve(_bins.size());
		y.reserve(_bins.size());
		for(std::map<int, Bin>::const_iterator i=_bins.begin();
				i!=_bins.end(); ++i)
		{
			double avg = i->second.average();
			x.push_back(binToValue(i->first));
			y.push_back(avg);
		}
	}
	
private:
	struct Bin
	{
		Bin() : sum(0.0), sumSq(0.0), weight(0.0) { }
		double sum, sumSq;
		double weight;
		double average() const { return sum / weight; }
		double stdErr() const
		{
			return sqrt(sumSq)/weight;
		}
		double sigmaSq() const
		{
			double sumMeanSquared = sum * sum / weight;
			return (sumSq - sumMeanSquared) / (weight-1.0);
		}
	};
	
	double _binsPerOrderOfMagnitude;
	std::map<int, Bin> _bins;
	
	int valueToBinCenter(double val) const
	{
		int centeredBin = valueCenterToBin(val);
		double leftBin = binToValue(centeredBin);
		double rightBin = binToValue(centeredBin+1);
		if(val - leftBin < rightBin - val)
			return centeredBin;
		else
			return centeredBin+1;
	}
	
	int valueCenterToBin(double val) const
	{
		if(val <= 0.0)
			return std::numeric_limits<int>::min();
		int bin = 0;
		while(val >= 10.0)
		{
			bin += _binsPerOrderOfMagnitude;
			val /= 10.0;
		}
		while(val < 1.0)
		{
			bin -= _binsPerOrderOfMagnitude;
			val *= 10.0;
		}
		double valStep = pow(10.0, 1.0/_binsPerOrderOfMagnitude);
		while(val >= valStep)
		{
			++bin;
			val /= valStep;
		}
		return bin;
	}
	
	double binToValue(int bin) const
	{
		int base = int(floor(double(bin)/_binsPerOrderOfMagnitude));
		bin -= base*_binsPerOrderOfMagnitude;
		double baseVal = exp10(base);
		double valStep = pow(10.0, 1.0/_binsPerOrderOfMagnitude);
		while(bin > 0)
		{
			baseVal *= valStep;
			--bin;
		}
		return baseVal;
	}
};

#endif
