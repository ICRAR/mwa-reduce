#include "rmsynthesis.h"
#include "banddata.h"

RMSynthesis::RMSynthesis(const SpectralEnergyDistribution& sed) :
	_sed(sed)
{
}

void RMSynthesis::Synthesize()
{
	SpectralEnergyDistribution::const_iterator t = _sed.begin();
	double minLambda = BandData::FrequencyToLambda(t->first);
	++t;
	double nextMinLambda = BandData::FrequencyToLambda(t->first);
	double maxLambda = BandData::FrequencyToLambda(_sed.rbegin()->first);
	size_t sampleCount = size_t(ceil(2.0 * maxLambda / (nextMinLambda - minLambda)));
	_fdf.assign(0.0, sampleCount);s
	for(SpectralEnergyDistribution::const_iterator i=_sed.begin(); i!=_sed.end(); ++i)
	{
		double
			lambda = BandData::FrequencyToLambda(i->first),
			q = i->second.FluxDensity(Polarization::StokesQ),
			u = i->second.FluxDensity(Polarization::StokesU);
		addSample(lambda, q, u);
	}
}

void RMSynthesis::addSample(double lambdaSq, double q, double u)
{
	for(size_t i=0; i!=_fdf.size(); ++i)
	{
		double x = double(i) * 2.0 * M_PI / _fdf.size();
		// (q + iu) exp ( i 2 pi x lambdasq ) + (q - iu) exp ( - i 2 pi x lambdasq )
		// = 2 real ( (q + iu) exp ( i 2 pi x lambdasq ) )
		// = 2 real ( (q + iu) [ cos ( 2 pi x lambdasq ) + i sin ( 2 pi x lambdasq ) ] )
		// = 2 [ q cos ( 2 pi x lambdasq ) - u sin ( 2 pi x lambdasq ) ]
		double s, c;
		sincos(2.0*M_PI*x*lambdaSq, &s, &c);
		_fdf[i] += 2.0 * (q * c - u * s);
	}
}
