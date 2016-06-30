#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <cmath>

class Universe
{
public:
	/**
	 * Mass density parameter
	 */
	static double OmegaM()
	{
		return 0.27;
	}
	
	/**
	 * Curvature
	 */
	static double OmegaK()
	{
		return 1.0 - OmegaM() - OmegaL();
	}
	
	/**
	 * Dark energy density parameter
	 * (sometimes OmegaV or lambda)
	 */
	static double OmegaL()
	{
		return 0.73;
	}
	
	
	static double H0()
	{
		return 0.72;
	}
	
	/**
	 * Hubble distance
	 */
	static double DH()
	{
		return C() / H0();
	}
	
	static double DHinMPC()
	{
		return 4228.0;
	}
	
	/**
	 * Speed of light
	 */
	static double C()
	{
		return 299792458.0;
	}
	
	/**
	 * Hubble parameter
	 */
	static double E(double z)
	{
		return sqrt(OmegaM()*(1.0+z)*(1.0+z)*(1.0+z) + OmegaK()*(1.0+z)*(1.0+z) + OmegaL());
	}
	
	static double OmegaR()
	{
		return 4.165e-5/(H0()*H0());
	}
	
	static double ComovingDistanceInMPCoverH(double z)
	{
		/*
		double d = 0.0;
		for(size_t i=0; i!=10001; ++i)
		{
			d += 1.0 / E(double(i)*z/10000.0);
		}
		// We integrate over 10001 samples :
		return d * DHinMPC() * z * H0() / 10001.0;*/
		
		const size_t n = 10000;
		double d = 0.0;
		double az = 1.0/(1.0+1.0*z);
		for(size_t i=0; i!=n; ++i)
		{
			double a = az+(1.0-az)*(double(i)+0.5)/double(n);
			double adot = sqrt(OmegaK()+(OmegaM()/a)+(OmegaR()/(a*a))+(OmegaL()*a*a));
			d += 1.0 / (a*adot);
		}
		
		d = (1.0-az)*d/n;
		return (C()/1000.0/(H0()*100.0)/*=h*/)*d*H0();
	}
	
	static double ComovingDistanceInMPC(double z)
	{
		/*
		double d = 0.0;
		for(size_t i=0; i!=10001; ++i)
		{
			d += 1.0 / E(double(i)*z/10000.0);
		}
		// We integrate over 10001 samples :
		return d * DHinMPC() * z * H0() / 10001.0;*/
		
		const size_t n = 10000;
		double d = 0.0;
		double az = 1.0/(1.0+1.0*z);
		for(size_t i=0; i!=n; ++i)
		{
			double a = az+(1.0-az)*(double(i)+0.5)/double(n);
			double adot = sqrt(OmegaK()+(OmegaM()/a)+(OmegaR()/(a*a))+(OmegaL()*a*a));
			d += 1.0 / (a*adot);
		}
		
		d = (1.0-az)*d/n;
		return (C()/1000.0/(H0()*100.0)/*=h*/)*d;
	}
	
	static double HIRestFrequencyInHz()
	{
		return 1420405752.0;
	}
	
	static double HIRedshift(double frequencyInHz)
	{
		return HIRestFrequencyInHz()/frequencyInHz - 1.0;
	}
	
	static double Boltzmann()
	{
		return 1.3806488e-23;
	}
};

#endif
