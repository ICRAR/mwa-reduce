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
			
			return SourceSDFWithSI<NumericType>(fluxA, freqA, fluxB, freqB).FluxAtFrequency(frequencyHz);
		}
			
		virtual NumericType IntegratedFlux(NumericType startFrequency, NumericType endFrequency) const
		{
			if(_fluxes.size() <= 2)
			{
				if(_fluxes.empty())
					return 0.0;
				else if(_fluxes.size()==1)
					return _fluxes.begin()->second * (endFrequency - startFrequency);
				else {
					NumericType
						freqA = _fluxes.begin()->first,
						fluxA = _fluxes.begin()->second,
						freqB = _fluxes.rbegin()->first,
						fluxB = _fluxes.rbegin()->second;
					
					return SourceSDFWithSI<NumericType>(fluxA, freqA, fluxB, freqB).IntegratedFlux(startFrequency, endFrequency);
				}
			}
				
			typename FluxMap::const_iterator iter = _fluxes.lower_bound(startFrequency);
			if(iter != _fluxes.begin()) --iter;
			
			NumericType leftFrequency = startFrequency;
			
			NumericType integratedSum = 0.0;
			while(iter != _fluxes.end() && iter->first < endFrequency)
			{
				typename FluxMap::const_iterator left = iter;
				typename FluxMap::const_iterator right = iter;
				++right;
				
				NumericType rightFrequency = right->first;
				// If this is past the sampled frequencies, extrapolate last 2 samples
				if(right == _fluxes.end()) {
					--right; --left;
					rightFrequency = endFrequency;
				}
				
				NumericType
					freqA = left->first,
					fluxA = left->second,
					freqB = right->first,
					fluxB = right->second;
				
				integratedSum += SourceSDFWithSI<NumericType>(fluxA, freqA, fluxB, freqB).IntegratedFlux(leftFrequency, rightFrequency);
				leftFrequency = rightFrequency;
				++iter;
			}
			return integratedSum;
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
			_fluxes.insert(std::pair<NumericType, NumericType>(frequency, flux));
		}
		
	private:
		FluxMap _fluxes;
};


#endif
