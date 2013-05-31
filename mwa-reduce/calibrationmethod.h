#ifndef CALIBRATION_METHOD_H
#define CALIBRATION_METHOD_H

#include "predicter.h"

#include <complex>
#include <cstddef>
#include <vector>

class CalibrationMethod
{
	private:
		template<typename T>
		class ValueArray
		{
			public:
				ValueArray(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
					_nBaselines(nAntenna * (nAntenna-1) / 2),
					_nChannels(nChannels),
					_nAntenna(nAntenna),
					_nTimesteps(nTimesteps),
					_values(nChannels * _nBaselines * nTimesteps)
				{
				}
				
				T *ValuePtr(size_t antenna1, size_t antenna2, size_t timeIndex)
				{
					return &_values[4 * _nChannels * (BaselineIndex(antenna1, antenna2) + timeIndex * _nBaselines)];
				}
			private:
				size_t BaselineIndex(size_t antenna1, size_t antenna2) const
				{
					return antenna1 * _nAntenna - (antenna1-1) * antenna1 / 2 + antenna2;
				}
				// Ordered in Polarization, Frequency, Antenna2, Antenna1 (i.e., Antenna1 is most significant, and a1 <= a2), Time
				size_t _nBaselines, _nChannels, _nAntenna, _nTimesteps;
				std::vector<T> _values;
		};
		
		typedef ValueArray<std::complex<double> > DataArray;
		typedef ValueArray<double> WeightArray;

	public:
		CalibrationMethod(size_t nChannels, size_t nAntenna, size_t nTimesteps);
		
		void AddData(const std::complex<float> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2, size_t timestep);
		
		void Execute(double precisionLimit);
		
	private:
		void calculateNextIter(size_t ant);
		
		static void multiplyWithInverse2x2(std::complex<double>* lhs, std::complex<double>* rhs);
		
		DataArray _data, _model;
		WeightArray _weights;
		size_t _nChannels, _nAntenna;
};

#endif
