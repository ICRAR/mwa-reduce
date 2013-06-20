#ifndef PREDICTER_H
#define PREDICTER_H

#include <complex>

class Predicter
{
	public:
		typedef double NumType;
		typedef std::complex<NumType> CNumType;
		
		Predicter(NumType phaseCentreRA, NumType phaseCentreDec, NumType startFrequency, NumType endFrequency, size_t channelCount) :
			_ra0(phaseCentreRA), _dec0(phaseCentreDec), _startFrequency(startFrequency), _endFrequency(endFrequency), _channelCount(channelCount)
		{
		}
		
		/**
		 * Initializes the l and m position(s) of the source.
		 */
		void Initialize(class ModelSource &model);
		void Initialize(class Model &model);
		
		CNumType Predict(const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
		CNumType Predict(const class Model &model, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
	private:
		struct SourceParameters
		{
			NumType l, m, lmsqrt, *brightness;
		};
		
		NumType _ra0, _dec0, _startFrequency, _endFrequency;
		size_t _channelCount;
};

#endif
