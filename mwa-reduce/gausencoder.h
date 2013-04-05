#ifndef GAUS_ENCODER_H
#define GAUS_ENCODER_H

#include <vector>
#include <cstring>
#include <algorithm>

template<typename ValueType=float>
class GausEncoder
{
	public:
		GausEncoder(size_t quantCount, ValueType stddev, bool gaussianMapping = true);
		
		typedef unsigned symbol_t;
		typedef ValueType value_t;
		
		symbol_t Encode(ValueType value) const
		{
			DictionaryItem item;
			item.value = value;
			typename std::vector<DictionaryItem>::const_iterator iter =
				std::lower_bound(_dictionary.begin(), _dictionary.end(), item);
			
			// if all keys are smaller, return largest symbol.
			if(iter == _dictionary.end())
				return _dictionary.rbegin()->symbol;
			// if smaller than first item return first item
			else if(iter == _dictionary.begin())
				return iter->symbol;
			
			// requested value is between two symbols: round to nearest.
			typename std::vector<DictionaryItem>::const_iterator prev = iter-1;
			if(value - prev->value < iter->value - value)
				return prev->symbol;
			else
				return iter->symbol;
		}
		
		ValueType Decode(symbol_t symbol) const
		{
			return _dictionary[symbol].value;
		}
		
	private:
		typedef long double num_t;
		
		static num_t cumulative(num_t x);
		static num_t invCumulative(num_t c, num_t err = num_t(1e-13));
		
		struct DictionaryItem
		{
			value_t value;
			symbol_t symbol;
			
			bool operator<(const DictionaryItem &other) const { return value < other.value; }
			bool operator<(ValueType other) const { return value < other; }
		};
		
		std::vector<DictionaryItem> _dictionary;
};

#endif
