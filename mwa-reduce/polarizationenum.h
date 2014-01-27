#ifndef POLARIZATION_ENUM_H
#define POLARIZATION_ENUM_H

#include <stdexcept>

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
		StokesV
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
};

typedef Polarization::PolarizationEnum PolarizationEnum;

#endif
