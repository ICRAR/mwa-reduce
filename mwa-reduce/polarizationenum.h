#ifndef POLARIZATION_ENUM_H
#define POLARIZATION_ENUM_H

#include <stdexcept>
#include <vector>

class Polarization
{
public:
	enum PolarizationEnum
	{
		XX,
		XY,
		YX,
		YY,
		StokesI,
		StokesQ,
		StokesU,
		StokesV,
		RR,
		RL,
		LR,
		LL
	};
	
	static size_t TypeToIndex(enum PolarizationEnum polarization, size_t polCountInSet)
	{
		switch(polCountInSet)
		{
			case 1:
				if(polarization != StokesI)
					throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
				else
					return 0;
			case 2:
			switch(polarization)
			{
				case XX: return 0;
				case YY: return 1;
				default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
			}
			case 4:
			switch(polarization)
			{
				case XX: return 0;
				case XY: return 1;
				case YX: return 2;
				case YY: return 3;
				default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
			}
			default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
		}
	}
	
	static size_t TypeTo4PolIndex(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			case XX: return 0;
			case XY: return 1;
			case YX: return 2;
			case YY: return 3;
			default: throw std::runtime_error("TypeTo4PolIndex(): can't convert given polarization to index");
		}
	}
	
	static std::string TypeToShortString(enum PolarizationEnum polarization)
	{
		switch(polarization)
		{
			case XX: return "XX";
			case XY: return "XY";
			case YX: return "YX";
			case YY: return "YY";
			case StokesI: return "I";
			case StokesQ: return "Q";
			case StokesU: return "U";
			case StokesV: return "V";
			case RR: return "RR";
			case RL: return "RL";
			case LR: return "LR";
			case LL: return "LL";
			default: return "";
		}
	}
	
	static std::vector<PolarizationEnum> ParseList(const std::string& listStr)
	{
		std::vector<PolarizationEnum> list;
		enum { StartSt, GotXSt, GotYSt, GotLSt, GotRSt, GotSeperatorSt } state = StartSt;
		for(std::string::const_iterator i=listStr.begin(); i!=listStr.end(); ++i)
		{
			char c = (*i>='a' && *i<='z') ? *i-('a'-'A') : *i;
			switch(c)
			{
				case 'X':
					if(state == StartSt || state == GotSeperatorSt) state=GotXSt;
					else {
						if(state==GotXSt) list.push_back(XX);
						else if(state==GotYSt) list.push_back(YX);
						else throw std::runtime_error("Invalid polarization list: parse error near 'X'");
						state=StartSt;
					}
					break;
				case 'Y':
					if(state == StartSt || state == GotSeperatorSt) state=GotYSt;
					else {
						if(state==GotXSt) list.push_back(XY);
						else if(state==GotYSt) list.push_back(YY);
						else throw std::runtime_error("Invalid polarization list: parse error near 'Y'");
						state=StartSt;
					}
					break;
				case 'R':
					if(state == StartSt || state == GotSeperatorSt) state=GotRSt;
					else {
						if(state==GotRSt) list.push_back(RR);
						else if(state==GotLSt) list.push_back(LR);
						else throw std::runtime_error("Invalid polarization list: parse error near 'R'");
						state=StartSt;
					}
					break;
				case 'L':
					if(state == StartSt || state == GotSeperatorSt) state=GotLSt;
					else {
						if(state==GotRSt) list.push_back(RL);
						else if(state==GotLSt) list.push_back(LL);
						else throw std::runtime_error("Invalid polarization list: parse error near 'L'");
						state=StartSt;
					}
					break;
				case 'I':
					if(state == StartSt || state == GotSeperatorSt) list.push_back(StokesI);
					else throw std::runtime_error("Invalid polarization list: parse error near 'I'");
					state = StartSt;
					break;
				case 'Q':
					if(state == StartSt || state == GotSeperatorSt) list.push_back(StokesQ);
					else throw std::runtime_error("Invalid polarization list: parse error near 'Q'");
					state = StartSt;
					break;
				case 'U':
					if(state == StartSt || state == GotSeperatorSt) list.push_back(StokesU);
					else throw std::runtime_error("Invalid polarization list: parse error near 'U'");
					state = StartSt;
					break;
				case 'V':
					if(state == StartSt || state == GotSeperatorSt) list.push_back(StokesV);
					else throw std::runtime_error("Invalid polarization list: parse error near 'V'");
					state = StartSt;
					break;
				case ',':
				case ' ':
				case '/':
					if(state == StartSt) state = GotSeperatorSt;
					else throw std::runtime_error("Invalid polarization list: parse error near seperator");
			}
		}
		if(state!=StartSt)
			throw std::runtime_error("Invalid polarization list: parse error near string end");
		return list;
	}
};

typedef Polarization::PolarizationEnum PolarizationEnum;

#endif
