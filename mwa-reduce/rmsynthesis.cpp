#include "rmsynthesis.h"
#include "banddata.h"

RMSynthesis::RMSynthesis(const SpectralEnergyDistribution& sed) :
	_sed(sed)
{
}

void RMSynthesis::Synthesize()
{
	SpectralEnergyDistribution::const_reverse_iterator t = _sed.rbegin();
	double minLSq = BandData::FrequencyToLambda(t->first) * BandData::FrequencyToLambda(t->first);
	++t;
	double nextMinLSq = BandData::FrequencyToLambda(t->first) * BandData::FrequencyToLambda(t->first);
	_maxLSq = BandData::FrequencyToLambda(_sed.begin()->first) * BandData::FrequencyToLambda(_sed.begin()->first);
	size_t sampleCount = size_t(ceil(2.0 * _maxLSq / (nextMinLSq - minLSq)));
	_fdf.assign(sampleCount, 0.0);
	for(SpectralEnergyDistribution::const_iterator i=_sed.begin(); i!=_sed.end(); ++i)
	{
		double
			lambda = BandData::FrequencyToLambda(i->first),
			q = i->second.FluxDensity(Polarization::StokesQ),
			u = i->second.FluxDensity(Polarization::StokesU);
		if(std::isfinite(q) && std::isfinite(u))
			addSample(lambda*lambda, q, u);
	}
}

void RMSynthesis::addSample(double lambdaSq, double q, double u)
{
	for(size_t i=0; i!=_fdf.size(); ++i)
	{
		double x = IndexToValue(i); //(double(i) - double(_fdf.size()/2))/100.0;
		// (q + iu) exp ( i 2 pi x lambdasq ) + (q - iu) exp ( - i 2 pi x lambdasq )
		// = 2 real ( (q + iu) exp ( i 2 pi x lambdasq ) )
		// = 2 real ( (q + iu) [ cos ( 2 pi x lambdasq ) + i sin ( 2 pi x lambdasq ) ] )
		// = 2 [ q cos ( 2 pi x lambdasq ) - u sin ( 2 pi x lambdasq ) ]
		double s, c;
		sincos(x*lambdaSq*2.0, &s, &c);
		_fdf[i] += 2.0 * (q * c + u * s);
	}
}
