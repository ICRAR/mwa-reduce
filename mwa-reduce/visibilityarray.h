#ifndef VISIBILITY_ARRAY_H
#define VISIBILITY_ARRAY_H

#include <stdexcept>
#include <sstream>

#include "uvector.h"

template<typename T, int NPol>
class VisibilityArray
{
public:
	VisibilityArray(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
		_nBaselines(nAntenna * (nAntenna-1) / 2),
		_nChannels(nChannels),
		_nAntenna(nAntenna),
		_nTimesteps(nTimesteps),
		_values(nChannels * _nBaselines * nTimesteps * NPol, T())
	{
	}
	
	T *ValuePtr(size_t antenna1, size_t antenna2, size_t timeIndex)
	{
		if(antenna1 == antenna2 || antenna1 >= _nAntenna || antenna2 >=_nAntenna || timeIndex >= _nTimesteps)
		{
			std::ostringstream errMsg;
			errMsg << "Call: ValuePtr(antenna1=" << antenna1 << ", antenna2=" << antenna2 << ", timeIndex=" << timeIndex << ") - ";
			std::string s = errMsg.str();
			if(antenna1 == antenna2)
				throw std::runtime_error(s+"Requested an auto-correlation");
			if(antenna1 >= _nAntenna)
				throw std::runtime_error(s+"Requested index value for antenna1 out of bounds");
			if(antenna2 >= _nAntenna)
				throw std::runtime_error(s+"Requested index value for antenna2 out of bounds");
			if(timeIndex >= _nTimesteps)
				throw std::runtime_error(s+"Requested index value for time index out of bounds");
		}
		if(antenna1 > antenna2)
			std::swap(antenna1, antenna2);
		return &_values[NPol * _nChannels * (BaselineIndex(antenna1, antenna2) + timeIndex * _nBaselines)];
	}
	void SetAll(const T& newValue)
	{
		_values.assign(_nChannels * _nBaselines * _nTimesteps * NPol, newValue);
	}
private:
	size_t BaselineIndex(size_t antenna1, size_t antenna2) const
	{
		return (antenna1*(2*_nAntenna - antenna1 - 3) + 2*antenna2 - 2)/2;
	}
	// Ordered in Polarization, Frequency, Antenna2, Antenna1 (i.e., Antenna1 is most significant, and a1 <= a2), Time
	size_t _nBaselines, _nPolarizations, _nChannels, _nAntenna, _nTimesteps;
	ao::uvector<T> _values;
};

#endif
