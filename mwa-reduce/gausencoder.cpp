#include "gausencoder.h"

#ifdef HAS_GSL
#include <gsl/gsl_sf_erf.h>
#endif

#include <stdexcept>
#include <limits>
#include <cmath>

template<typename ValueType>
inline typename GausEncoder<ValueType>::num_t GausEncoder<ValueType>::cumulative(num_t x)
{
#ifdef HAS_GSL
	return num_t(0.5) + num_t(0.5) * gsl_sf_erf(x/num_t(M_SQRT2l));
#else
	throw std::runtime_error("GausEncoder was compiled without GSL: can't calculate cumulative distribution without GSL");
#endif
}

template<typename ValueType>
typename GausEncoder<ValueType>::num_t GausEncoder<ValueType>::invCumulative(num_t c, num_t err)
{
	if(c < 0.5) return(-invCumulative(1.0 - c, err));
	else if(c == 0.5) return 0.0;
	else if(c == 1.0) return std::numeric_limits<num_t>::infinity();
	else if(c > 1.0) return std::numeric_limits<num_t>::quiet_NaN();
	
	num_t x = 1.0;
	num_t fx = cumulative(x);
	num_t xLow, xHi;
	if(fx < c)
	{
		do {
			x *= 2.0;
			fx = cumulative(x);
		} while(fx < c);
		xLow = x * 0.5;
		xHi = x;
	} else {
		xLow = 0.0;
		xHi = 1.0;
	}
	num_t error = xHi;
	int notConverging = 0;
	do {
		x = (xLow + xHi) * 0.5;
		fx = cumulative(x);
		if(fx > c)
			xHi = x;
		else
			xLow = x;
		num_t currErr = std::fabs(fx - c);
		if(currErr >= error)
		{
			++notConverging;
			// not converging anymore; stop.
			if(notConverging > 10)
				return x;
		} else
			notConverging = 0;
		error = currErr;
	} while(error > err);
	return x;
}

template<typename ValueType>
GausEncoder<ValueType>::GausEncoder(size_t quantCount, ValueType stddev, bool gaussianMapping) :
	_dictionary(quantCount)
{
	if(gaussianMapping)
	{
		for(size_t i=0; i!=quantCount; ++i)
		{
			DictionaryItem &item = _dictionary[i];
			num_t val = ((num_t) i + num_t(0.5)) / (num_t) (quantCount);
			item.value = stddev * invCumulative(val);
			item.symbol = i;
		}
	} else {
		for(size_t i=0; i!=quantCount; ++i)
		{
			DictionaryItem &item = _dictionary[i];
			num_t val = -1.0 + 2.0*((num_t) i + num_t(0.5)) / (num_t) (quantCount);
			item.value = stddev * val;
			item.symbol = i;
		}
	}
}

template class GausEncoder<float>;

