#ifndef ANGLE_H
#define ANGLE_H

#include <string>
#include <stdexcept>
#include <cmath>

#include <boost/algorithm/string/case_conv.hpp>

class Angle
{
public:
	/**
	 * Parse the string as an angle, possibly with unit specification, and return in radians.
	 * @return The angle
	 */
	static double Parse(const std::string& s, const std::string& valueDescription);
	
private:
	static size_t findNumberEnd(const std::string& s);
	static bool isDigit(const char c) { return c>='0' && c<='9'; }
	static bool isWhitespace(const char c) { return c==' ' || c=='\t'; }
};

inline double Angle::Parse(const std::string& s, const std::string& valueDescription)
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

	// In degrees? (default)
	if(unitStr.empty() || unitStr=="deg" || unitStr=="degrees")
		return val * M_PI/180.0;
	
	// In arcmin?
	else if(unitStr.empty() || unitStr=="amin" || unitStr=="arcmin" || unitStr=="\'")
		return val * M_PI/(180.0*60.0);
	
	// In arcsec?
	else if(unitStr.empty() || unitStr=="asec" || unitStr=="arcsec" || unitStr=="\'\'")
		return val * M_PI/(180.0*60.0*60.0);
	
	// In radians
	else if(unitStr.empty() || unitStr=="rad" || unitStr=="radians")
		return val * M_PI/(180.0*60.0*60.0);
	
	else
		throw std::runtime_error("Invalid unit specification in angle given for " + valueDescription);
}

inline size_t Angle::findNumberEnd(const std::string& s)
{
	const char* c = s.c_str();
	size_t pos = 0;
	while(isWhitespace(c[pos]))
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
