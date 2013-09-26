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
	typedef ao::uvector<std::pair<double, double>> WeightMap;
public:
	UvwDistribution(size_t binCount = 100) : _binCount(binCount)
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
		_maxBaseline = sqrt(maxDistSq);
		_minBaseline = sqrt(minDistSq);
		_maxUVW = _maxBaseline / bandData.SmallestWavelength();
		_minUVW = _minBaseline / bandData.LongestWavelength();
		const double
			uvwRange = _maxUVW - _minUVW,
			baselineRange = _maxBaseline - _minBaseline;
			
		ao::uvector<size_t> uvwHist(_binCount, 0), baselineLengthHist(_binCount, 0);
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
			
			casa::Vector<double> uvwVec = uvwColumn(row);
			const double distInM = sqrt(distSq(uvwVec));
			if(!skip)
			{
				for(size_t ch=0; ch!=bandData.ChannelCount(); ++ch)
				{
					double uvwDist = distInM / bandData.ChannelWavelength(ch);
					double bin = round((uvwDist - _minUVW) * double(_binCount-1) / uvwRange);
					size_t binIndex = size_t(bin);
					if(binIndex < _binCount)
					{
						uvwHist[size_t(bin)]++;
					}
				}
			}
			
			double bbin = round((distInM - _minBaseline) * double(_binCount-1) / baselineRange);
			size_t bbinIndex = size_t(bbin);
			if(bbinIndex < _binCount)
			{
				baselineLengthHist[size_t(bbin)]++;
			}
		}
		
		size_t maxCount = *std::max_element(uvwHist.begin(), uvwHist.end());
		double weightFactor = maxCount;
		_normalizationFactor = 1.0 / maxCount;
		
		_weightMap.resize(uvwHist.size());
		_baselineLengthMap.resize(uvwHist.size());
		for(size_t i=0; i!=uvwHist.size(); ++i)
		{
			double binCentre = (double(i) * uvwRange / double(_binCount) + _minUVW);
			double weight = weightFactor / double(uvwHist[i]);
			_weightMap[i] = std::make_pair(binCentre, weight);
			
			double bbinCentre = (double(i) * baselineRange / double(_binCount) + _minBaseline);
			_baselineLengthMap[i] = std::make_pair(bbinCentre, baselineLengthHist[i]*2);
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
	
	double CumulativeCount(double baselineDist) const
	{
		double count = 0.0;
		const WeightMap::const_iterator upper = std::lower_bound(_baselineLengthMap.begin(), _baselineLengthMap.end(), std::make_pair(baselineDist, 0.0));
		for(WeightMap::const_iterator i=upper; i!=_baselineLengthMap.end(); ++i)
		{
			count += i->second;
		}
		return count;
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
	double MinBaseline() const { return _minBaseline; }
	double MaxBaseline() const { return _maxBaseline; }
	size_t BinCount() const { return _binCount; }
private:
	const size_t _binCount;
	double _minUVW, _maxUVW;
	double _minBaseline, _maxBaseline;
	double _normalizationFactor;
	double _plFactor, _plExp;
	ao::uvector<std::pair<double, double>> _weightMap;
	ao::uvector<std::pair<double, double>> _baselineLengthMap;
	
	double distSq(casa::Vector<double>& uvwVec)
	{
		double u = uvwVec[0], v = uvwVec[1], w = uvwVec[2];
		return u*u + v*v + w*w;
	}
};

#endif
