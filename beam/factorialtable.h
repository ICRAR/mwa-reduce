#ifndef FACTORIAL_TABLE_H
#define FACTORIAL_TABLE_H

#include <boost/math/special_functions/factorials.hpp>

#include <aocommon/uvector.h>

class FactorialTable
{
public:
	FactorialTable(size_t nPrecalculated) : _table(nPrecalculated)
	{
		for(unsigned i=0;i<=nPrecalculated;i++)
			_table[i] = boost::math::factorial<double>(i);
	}
	
	double operator()(unsigned n) const
	{
		if(n >= _table.size())
		{
			return boost::math::factorial<double>(n);
		}
		
		if(n < _table.size())
			return _table[n];
		else
			return boost::math::factorial<double>(n);
	}	
	
private:
	aocommon::UVector<double> _table;
};

#endif
