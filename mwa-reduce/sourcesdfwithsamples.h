#include "sourcesdf.h"

#ifndef SOURCE_SDF_WITH_SAMPLES
#define SOURCE_SDF_WITH_SAMPLES

#include <map>

template<typename NumericType=long double>
class SourceSDFWithSamples : public SourceSDF<NumericType>
{
	private:
		typedef std::map<NumericType, NumericType> FluxMap;
	public:
		SourceSDFWithSamples() { }
		
		SourceSDFWithSamples(const SourceSDFWithSamples<NumericType> &source)
		: _fluxes(source._fluxes)
		{ }
		
		virtual NumericType FluxAtFrequency(NumericType frequencyHz) const
		{
			if(_fluxes.size() <= 1)
			{
				if(_fluxes.empty()) return 0.0;
				else return _fluxes.begin()->second;
			}
			
			// 'right' will be first item which frequency >= frequencyHz
			typename FluxMap::const_iterator right = _fluxes.lower_bound(frequencyHz);
			if(right->first == frequencyHz)
				return right->second;
			
			typename FluxMap::const_iterator left;
			
			// If the requested frequency is outside the range, we extrapolate the SI of the full range
			if(right == _fluxes.begin() || right == _fluxes.end())
			{
				left = _fluxes.begin();
				right = _fluxes.end();
				--right;
			} else {
				// Requested frequency is within range (no extrapolation required)
				left = right;
				--left;
			}
			
			NumericType
				freqA = left->first,
				fluxA = left->second,
				freqB = right->first,
				fluxB = right->second;
			
			return SourceSDFWithSI<NumericType>::FluxAtFrequency(fluxA, freqA, fluxB, freqB, frequencyHz);
		}
			
		virtual NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const
		{
			if(startFrequency == endFrequency)
				return FluxAtFrequency(startFrequency);
			
			typename FluxMap::const_iterator iter = _fluxes.lower_bound(startFrequency);
			
			/** Handle special cases */
			if(_fluxes.size() <= 2)
			{
				if(_fluxes.empty())
					return 0.0;
				else if(_fluxes.size()==1)
					return _fluxes.begin()->second;
				else if(_fluxes.size()==2) {
					NumericType
						freqA = _fluxes.begin()->first,
						fluxA = _fluxes.begin()->second,
						freqB = _fluxes.rbegin()->first,
						fluxB = _fluxes.rbegin()->second;
					
					return SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
				}
			}
			if(iter == _fluxes.end()) { // all keys are lower, so take entire range
				NumericType
					freqA = _fluxes.begin()->first,
					fluxA = _fluxes.begin()->second,
					freqB = _fluxes.rbegin()->first,
					fluxB = _fluxes.rbegin()->second;
				return SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
			
			if(iter != _fluxes.begin()) --iter;
			
			if(iter->first >= endFrequency) {
				// all keys are outside range, higher than range
				NumericType
					freqA = _fluxes.begin()->first,
					fluxA = _fluxes.begin()->second,
					freqB = _fluxes.rbegin()->first,
					fluxB = _fluxes.rbegin()->second;
				return SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
			
			NumericType integratedSum = 0.0;
			NumericType leftFrequency = startFrequency;
			if(leftFrequency < iter->first)
			{
				// requested frequency is below first item; extrapolate
				NumericType
					freqA = _fluxes.begin()->first,
					fluxA = _fluxes.begin()->second,
					freqB = _fluxes.rbegin()->first,
					fluxB = _fluxes.rbegin()->second;
				NumericType sumTerm = SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, iter->first);
				integratedSum += sumTerm * (iter->first - startFrequency);
				leftFrequency = iter->first;
			}
				
			while(iter != _fluxes.end() && iter->first < endFrequency)
			{
				typename FluxMap::const_iterator left = iter;
				typename FluxMap::const_iterator right = iter;
				++right;
				
				NumericType rightFrequency;
				
				// If this is past the sampled frequencies, extrapolate full range
				if(right == _fluxes.end()) {
					left = _fluxes.begin();
					right = _fluxes.end();
					--right;
					rightFrequency = endFrequency;
				} else {
					rightFrequency = right->first;
					if(rightFrequency > endFrequency)
						rightFrequency = endFrequency;
				}
				
				NumericType
					freqA = left->first,
					fluxA = left->second,
					freqB = right->first,
					fluxB = right->second;
				
				if(leftFrequency < rightFrequency)
				{
					NumericType sumTerm = SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, leftFrequency, rightFrequency);
					if(!std::isfinite(sumTerm))
					{
						std::cerr << "Warning: integrating flux between " << leftFrequency << " and " << rightFrequency << " with fluxes " << fluxA << '@' << freqA << ',' << fluxB << '@' << freqB << " gave non-finite result\n";
					}
					
					integratedSum += sumTerm * (rightFrequency - leftFrequency);
				}
				leftFrequency = rightFrequency;
				++iter;
			}
			return integratedSum / (endFrequency - startFrequency);
		}
			
		virtual SourceSDF<NumericType> *Copy() const
		{
			SourceSDFWithSamples<NumericType> *sdf = new SourceSDFWithSamples<NumericType>();
			sdf->_fluxes = _fluxes;
			return sdf;
		}
		
		virtual std::string ToString() const
		{
			std::ostringstream s;
			s << "sampled "
			  << _fluxes.size();
			for(typename FluxMap::const_iterator i=_fluxes.begin(); i!=_fluxes.end(); ++i)
			  s << ' ' << (i->second) << ' ' << (i->first/1000000.0);
			return s.str();
		}
		
		void AddSample(NumericType flux, NumericType frequency)
		{
			if(std::isfinite(flux))
				_fluxes.insert(std::pair<NumericType, NumericType>(frequency, flux));
			//else
			//	std::cerr << "Warning: ignoring non-finite result\n";
		}
		
		typedef typename std::map<NumericType, NumericType>::const_iterator const_iterator;
		typedef typename std::map<NumericType, NumericType>::iterator iterator;
		
		const_iterator begin() const
		{
			return _fluxes.begin();
		}
		iterator begin()
		{
			return _fluxes.begin();
		}
		const_iterator end() const
		{
			return _fluxes.end();
		}
		iterator end()
		{
			return _fluxes.end();
		}
	private:
		FluxMap _fluxes;
};


#endif
