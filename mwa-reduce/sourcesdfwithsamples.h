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
			
			typename FluxMap::const_iterator after = _fluxes.lower_bound(frequencyHz);
			if(after->first == frequencyHz)
				return after->second;
			
			// If this frequency is outside the range, we extrapolate the nearest SI
			if(after == _fluxes.begin())
			{
				++after;
			} else if(after == _fluxes.end())
			{
				--after;
			}
			
			typename FluxMap::const_iterator before = after;
			--before;
			
			NumericType
				freqA = before->first,
				fluxA = before->second,
				freqB = after->first,
				fluxB = after->second;
			
			return SourceSDFWithSI<NumericType>::FluxAtFrequency(fluxA, freqA, fluxB, freqB, frequencyHz);
		}
			
		virtual NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const
		{
			if(startFrequency == endFrequency)
				return FluxAtFrequency(startFrequency);
			
			typename FluxMap::const_iterator iter = _fluxes.lower_bound(startFrequency);
			
			/** Handle special cases */
			if(_fluxes.size() <= 2 || iter == _fluxes.end())
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
				} else { // all keys are lower, so take last two
					typename FluxMap::const_reverse_iterator end = _fluxes.rbegin();
					NumericType
						freqB = end->first,
						fluxB = end->second;
						++end;
					NumericType
						freqA = end->first,
						fluxA = end->second;
					return SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
				}
			}
			
			if(iter != _fluxes.begin()) --iter;
			
			if(iter->first >= endFrequency) {
				// all keys are outside range, higher than range
				typename FluxMap::const_iterator begin = _fluxes.begin();
				NumericType
					freqA = begin->first,
					fluxA = begin->second;
					++begin;
				NumericType
					freqB = begin->first,
					fluxB = begin->second;
				return SourceSDFWithSI<NumericType>::IntegratedFlux(fluxA, freqA, fluxB, freqB, startFrequency, endFrequency);
			}
				
			NumericType leftFrequency = startFrequency;
			
			NumericType integratedSum = 0.0;
			while(iter != _fluxes.end() && iter->first < endFrequency)
			{
				typename FluxMap::const_iterator left = iter;
				typename FluxMap::const_iterator right = iter;
				++right;
				
				NumericType rightFrequency;
				// If this is past the sampled frequencies, extrapolate last 2 samples
				if(right == _fluxes.end()) {
					--right; --left;
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
		
		typename std::map<NumericType, NumericType>::const_iterator begin() const
		{
			return _fluxes.begin();
		}
		
		typename std::map<NumericType, NumericType>::const_iterator end() const
		{
			return _fluxes.end();
		}
	private:
		FluxMap _fluxes;
};


#endif
