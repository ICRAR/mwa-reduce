#ifndef ANGLE_H
#define ANGLE_H

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include <boost/algorithm/string/case_conv.hpp>

class Angle
{
public:
	enum Unit { Radians, Degrees, Arcminutes, Arcseconds, Milliarcseconds };
	/**
	 * Parse the string as an angle, possibly with unit specification, and return in radians.
	 * @return The angle
	 */
	static double Parse(const std::string& s, const std::string& valueDescription, Unit defaultUnit);
	
	static std::string ToNiceString(double angleRad);
	
private:
	static size_t findNumberEnd(const std::string& s);
	static bool isDigit(const char c) { return c>='0' && c<='9'; }
	static bool isWhitespace(const char c) { return c==' ' || c=='\t'; }
};

inline std::string Angle::ToNiceString(double angleRad)
{
	std::ostringstream str;
	double degAngle = angleRad * 180.0 / M_PI;
	if(degAngle < 0.0)
	{
		str << "-";
		degAngle = -degAngle;
	}
	if(degAngle >= 2.0)
	{
		str << round(degAngle*100.0)/100.0 << " deg";
	}
	else {
		double minAngle = angleRad * 180.0 * 60.0 / M_PI;
		if(minAngle >= 2.0)
		{
			str << round(minAngle*100.0)/100.0 << "'";
		}
		else {
			double secAngle = angleRad * 180.0 * 60.0 * 60.0 / M_PI;
			if(secAngle >= 1.0)
			{
				str << round(secAngle*100.0)/100.0 << "''";
			}
			else {
				str << round (secAngle*100.0*1000.0)/100.0 << " masec";
			}
		}
	}
	return str.str();
}

inline double Angle::Parse(const std::string& s, const std::string& valueDescription, Unit defaultUnit)
{
  size_t end = findNumberEnd(s);
	if(end == 0)
		throw std::runtime_error("Error parsing " + valueDescription);
	std::string number = s.substr(0, end);
	double val = atof(number.c_str());
	// Skip whitespace after number
	const char *c = s.c_str();
	while(isWhitespace(c[end]))
		++end;
	std::string unitStr = std::string(&c[end]);
	boost::to_lower(unitStr);

	// Unit string empty? Than use default unit.
	if(unitStr.empty())
	{
		switch(defaultUnit)
		{
			case Radians:
				return val;
			case Degrees:
				return val * M_PI/180.0;
			case Arcminutes:
				return val * M_PI/(180.0*60.0);
			case Arcseconds:
				return val * M_PI/(180.0*60.0*60.0);
			case Milliarcseconds:
				return val * M_PI/(180.0*60.0*60.0*1000.0);
		}
	}
	
	// In degrees?
	else if(unitStr=="deg" || unitStr=="degrees")
		return val * M_PI/180.0;
	
	// In arcmin?
	else if(unitStr=="amin" || unitStr=="arcmin" || unitStr=="\'")
		return val * M_PI/(180.0*60.0);
	
	// In arcsec?
	else if(unitStr=="asec" || unitStr=="arcsec" || unitStr=="\'\'")
		return val * M_PI/(180.0*60.0*60.0);
	
	// In marcsec?
	else if(unitStr=="mas" || unitStr=="masec" || unitStr=="marcsec")
		return val * M_PI/(180.0*60.0*60.0*1000.0);
	
	// In radians
	else if(unitStr=="rad" || unitStr=="radians")
		return val;
	
	throw std::runtime_error("Invalid unit specification in angle given for " + valueDescription);
}

inline size_t Angle::findNumberEnd(const std::string& s)
{
	const char* c = s.c_str();
	size_t pos = 0;
	while(isWhitespace(c[pos]))
		++pos;
	if(c[pos] == '-')
		++pos;
	while(isDigit(c[pos]))
		++pos;
	if(c[pos]=='.')
		++pos;
	while(isDigit(c[pos]))
		++pos;
	if(c[pos]=='e' || c[pos]=='E')
	{
		++pos;
		if(c[pos]=='-' || c[pos]=='+')
			++pos;
		while(isDigit(c[pos]))
			++pos;
	}
	return pos;
}

#endif
