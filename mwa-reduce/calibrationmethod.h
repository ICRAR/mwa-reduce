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
				ValueArray(size_t nChannels, size_t nAntenna) :
					_values(nChannels * nAntenna * (nAntenna-1)),
					_nChannels(nChannels),
					_nAntenna(nAntenna)
				{
				}
				
				T *ValuePtr(size_t antenna1, size_t antenna2)
				{
					return &_values[4 * _nChannels * BaselineIndex(antenna1, antenna2)];
				}
			private:
				size_t BaselineIndex(size_t antenna1, size_t antenna2) const
				{
					return antenna1 * _nAntenna - (antenna1-1) * antenna1 / 2 + antenna2;
				}
				// Ordered in Polarization, Frequency, Antenna2, Antenna1 (i.e., Antenna1 is most significant, and a1 <= a2)
				std::vector<T> _values;
				size_t _nChannels, _nAntenna;
		};
		
		typedef ValueArray<std::complex<double> > DataArray;
		typedef ValueArray<double> WeightArray;

	public:
		CalibrationMethod(size_t nChannels, size_t nAntenna);
		
		void AddData(const std::complex<float> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2);
		
		void Execute(double precisionLimit);
		
	private:
		static void multiplyWithInverse2x2(std::complex<double>* lhs, std::complex<double>* rhs);
		
		DataArray _data;
		WeightArray _weights;
		size_t _nChannels, _nAntenna;
};

#endif
