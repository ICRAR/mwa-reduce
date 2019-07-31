#ifndef WINDOW_FUNCTION_H
#define WINDOW_FUNCTION_H

#include <cmath>
#include <cstring>

class WindowFunction
{
public:
	enum Type { Rectangular, BlackmanNutall, BlackmanHarris, Hann, Tukey };
	
	static Type GetType(const std::string& typeStr)
	{
		if(typeStr == "hann")
			return Hann;
		else if(typeStr == "blackman-nutall")
			return BlackmanNutall;
		else if(typeStr == "blackman-harris")
			return BlackmanHarris;
		else if(typeStr == "rectangular")
			return Rectangular;
		else if(typeStr == "tukey")
			return Tukey;
		else
			throw std::runtime_error("The window function name is not a valid. Valid windows are: rectangular, hann, blackman-harris or blackman-nutall");
	}

	static double Evaluate(Type windowFunc, size_t n, size_t i)
	{
		switch(windowFunc)
		{
			case Rectangular: return 1.0;
			case BlackmanNutall: return EvaluateBlackmanNutall(n, i);
			case BlackmanHarris: return EvaluateBlackmanHarris(n, i);
			case Tukey: throw std::runtime_error("Tukey window requires parameter");
			case Hann: return EvaluateHann(n, i);
		}
		return 0.0;
	}
	
	static double EvaluateBlackmanNutall(size_t n, size_t i)
	{
		const static double
			a0=0.3635819, a1=0.4891775, a2=0.1365995, a3=0.0106411;
		const double
			id = double(i) * 2.0 * M_PI, nd = int(n);
		return
			a0 -
			a1 * std::cos( (   id   )/nd ) +
			a2 * std::cos( (2.0 * id)/nd ) -
			a3 * std::cos( (3.0 * id)/nd );
	}
	
	static double EvaluateBlackmanHarris(size_t n, size_t i)
	{
		const static double
			a0=0.35875, a1=0.48829, a2=0.14128, a3=0.01168;
		const double
			id = double(i) * 2.0 * M_PI, nd = int(n);
		return
			a0 -
			a1 * std::cos( (   id   )/nd ) +
			a2 * std::cos( (2.0 * id)/nd ) -
			a3 * std::cos( (3.0 * id)/nd );
	}
	
	static double EvaluateHann(size_t n, size_t i)
	{
		double s = std::sin(M_PI * double(i)/(double(n)));
		return s*s;
	}
};

#endif
