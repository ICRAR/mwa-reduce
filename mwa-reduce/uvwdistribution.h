#ifndef UVW_DISTRIBUTION_H
#define UVW_DISTRIBUTION_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>

#include "uvector.h"
#include "banddata.h"
#include "nlplfitter.h"

class UvwDistribution
{
private:
	typedef uvector<std::pair<double, double>> WeightMap;
public:
	UvwDistribution() : _binCount(100)
	{
	}
	
	void Calculate(casa::MeasurementSet& ms)
	{
		bool skipAutos = true;
		
		casa::ROScalarColumn<int> antenna1Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> antenna2Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
		casa::ROArrayColumn<double> uvwColumn(ms, casa::MS::columnName(casa::MSMainEnums::UVW));
		BandData bandData(ms.spectralWindow());
		
		bool isFirstValue = true;
		double minDistSq = 0.0, maxDistSq = 0.0;
		const double startTime = timeColumn(0);
		for(size_t row=0; row!=ms.nrow(); ++row)
		{
			if(startTime != timeColumn(row))
				break;
			
			bool skip = false;
			if(skipAutos)
			{
				const int
					a1 = antenna1Column(row),
					a2 = antenna2Column(row);
				if(a1 == a2)
					skip = true;
			}
			
			if(!skip)
			{
				casa::Vector<double> uvwVec = uvwColumn(row);
				const double dSq = distSq(uvwVec);
				if(isFirstValue)
				{
					minDistSq = dSq;
					maxDistSq = dSq;
					isFirstValue = false;
				}
				else {
					if(dSq > maxDistSq)
						maxDistSq = dSq;
					if(dSq < minDistSq)
						minDistSq = dSq;
				}
			}
		}
		_maxUVW = sqrt(maxDistSq) / bandData.SmallestWavelength();
		_minUVW = sqrt(minDistSq) / bandData.LongestWavelength();
		const double uvwRange = _maxUVW - _minUVW;
			
		uvector<size_t> hist(_binCount, 0);
		for(size_t row=0; row!=ms.nrow(); ++row)
		{
			if(startTime != timeColumn(row))
				break;
			
			bool skip = false;
			if(skipAutos)
			{
				const int
					a1 = antenna1Column(row),
					a2 = antenna2Column(row);
				if(a1 == a2)
					skip = true;
			}
			
			if(!skip)
			{
				casa::Vector<double> uvwVec = uvwColumn(row);
				const double distInM = sqrt(distSq(uvwVec));
				for(size_t ch=0; ch!=bandData.ChannelCount(); ++ch)
				{
					double uvwDist = distInM / bandData.ChannelWavelength(ch);
					double bin = round((uvwDist - _minUVW) * double(_binCount-1) / uvwRange);
					size_t binIndex = size_t(bin);
					if(binIndex < _binCount)
					{
						hist[size_t(bin)]++;
					}
				}
			}
		}
		
		size_t maxCount = *std::max_element(hist.begin(), hist.end());
		double weightFactor = maxCount;
		
		_weightMap.resize(hist.size());
		for(size_t i=0; i!=hist.size(); ++i)
		{
			double binCentre = (double(i) * uvwRange / double(_binCount) + _minUVW);
			double weight = weightFactor / double(hist[i]);
			_weightMap[i] = std::make_pair(binCentre, weight);
		}
		
		FitPowerlaw(_plExp, _plFactor);
	}
	
	double CountWithInterpolation(double uvwDistance) const
	{
		return 1.0 / WeightWithInterpolation(uvwDistance);
	}
	
	double WeightWithInterpolation(double uvwDistance) const
	{
		const WeightMap::const_iterator upper = std::lower_bound(_weightMap.begin(), _weightMap.end(), std::make_pair(uvwDistance, 0.0));
		if(upper == _weightMap.end())
			return _weightMap.back().second;
		if(upper == _weightMap.begin())
			return _weightMap.front().second;
		const WeightMap::const_iterator lower = upper-1;
		const double leftBin = lower->first, rightBin = upper->first;
		// Interpolate linearly between the two bins
		const double ratio = (uvwDistance - leftBin) / (rightBin - leftBin);
		return (1.0 - ratio) * lower->second + ratio * upper->second;
	}
	
	double WeightFromFit(double uvwDistance) const
	{
		return exp(_plFactor * pow(uvwDistance, _plExp));
	}
	
	void FitPowerlaw(double& exponent, double& factor) const
	{
		NonLinearPowerLawFitter fitter;
		for(WeightMap::const_iterator i=_weightMap.begin(); i!=_weightMap.end(); ++i)
		{
			fitter.AddDataPoint(i->first, log(i->second));
		}
		fitter.FastFit(exponent, factor);
	}
	
	double MinUvw() const { return _minUVW; }
	double MaxUvw() const { return _maxUVW; }
	size_t BinCount() const { return _binCount; }
private:
	size_t _binCount;
	double _minUVW, _maxUVW;
	double _plFactor, _plExp;
	uvector<std::pair<double, double>> _weightMap;
	
	double distSq(casa::Vector<double>& uvwVec)
	{
		double u = uvwVec[0], v = uvwVec[1], w = uvwVec[2];
		return u*u + v*v + w*w;
	}
};

#endif
