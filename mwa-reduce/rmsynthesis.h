#ifndef RM_SYNTHESIS_H
#define RM_SYNTHESIS_H

#include "spectralenergydistribution.h"
#include "uvector.h"

class RMSynthesis
{
public:
	RMSynthesis(const SpectralEnergyDistribution& sed);
	
	void Synthesize();
	
private:
	void addSample(double lambdaSq, double q, double u);
	
	SpectralEnergyDistribution _sed;
	ao::uvector<double> _fdf;
};

#endif
