#ifndef RM_SYNTHESIS_H
#define RM_SYNTHESIS_H

#include "model/measuredsed.h"

#include "uvector.h"

class RMSynthesis
{
public:
	explicit RMSynthesis(const MeasuredSED& sed);
	
	void Synthesize();
	
	const ao::uvector<double> &FDF() const { return _fdf; }
	
	double IndexToValue(size_t i) const { return (double(i) - double(_fdf.size()/2)) *0.5; } // TODO
	
	void PeakWithSmallRM(double maxRM, double& rmPeak, double& peakValue)
	{
		peakValue = 0.0;
		rmPeak = 0.0;
		for(size_t i=0; i!=_fdf.size(); ++i)
		{
			double rm = IndexToValue(i);
			if(rm <= maxRM)
			{
				double val = std::fabs(_fdf[i]);
				if(val > peakValue)
				{
					peakValue = val;
					rmPeak = rm;
				}
			}
		}
	}
	
	void PeakWithNonZeroRM(double minRM, double& rmPeak, double& peakValue)
	{
		peakValue = 0.0;
		rmPeak = 0.0;
		for(size_t i=0; i!=_fdf.size(); ++i)
		{
			double rm = IndexToValue(i);
			if(std::fabs(rm) >= minRM)
			{
				double val = std::fabs(_fdf[i]);
				if(val > peakValue)
				{
					peakValue = val;
					rmPeak = rm;
				}
			}
		}
	}
private:
	void addSample(double lambdaSq, double q, double u);
	
	MeasuredSED _sed;
	ao::uvector<double> _fdf;
	double _maxLSq;
};

#endif
