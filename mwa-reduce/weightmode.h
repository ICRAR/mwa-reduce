#ifndef WEIGHTMODE_H
#define WEIGHTMODE_H

#include <string>
#include <sstream>

class WeightMode
{
public:
	enum WeightingEnum { NaturalWeighted, DistanceWeighted, UniformWeighted, BriggsWeighted };
		
	WeightMode(WeightingEnum mode) : _mode(mode), _briggsRobustness(0.0)
	{ }
	
	WeightMode(const WeightMode& source) : _mode(source._mode), _briggsRobustness(source._briggsRobustness)
	{ }
	
	WeightMode& operator=(const WeightMode& source)
	{
		_mode = source._mode;
		_briggsRobustness = source._briggsRobustness;
		return *this;
	}
	
	bool operator==(const WeightMode& rhs)
	{
		if(_mode != rhs._mode)
			return false;
		else if(_mode == BriggsWeighted)
			return _briggsRobustness == rhs._briggsRobustness;
		else
			return true;
	}
	
	static WeightMode Briggs(double briggsRobustness)
	{
		WeightMode m(BriggsWeighted);
		m._briggsRobustness = briggsRobustness;
		return m;
	}
	
	WeightingEnum Mode() const { return _mode; }
	bool IsNatural() const { return _mode == NaturalWeighted; }
	bool IsDistance() const { return _mode == DistanceWeighted; }
	bool IsUniform() const { return _mode == UniformWeighted; }
	bool IsBriggs() const { return _mode == BriggsWeighted; }
	double BriggsRobustness() const { return _briggsRobustness; }
	bool RequiresGridding() const { return IsUniform() || IsBriggs(); }
	
	std::string ToString() const 
	{
		switch(_mode)
		{
			case UniformWeighted: return "uniform";
			case DistanceWeighted: return "distance";
			case NaturalWeighted: return "natural";
			case BriggsWeighted:
			{
				std::ostringstream s;
				s << "Briggs'(" << _briggsRobustness << ")";
				return s.str();
			}
		}
	}
private:
	enum WeightingEnum _mode;
	double _briggsRobustness;
};

#endif
