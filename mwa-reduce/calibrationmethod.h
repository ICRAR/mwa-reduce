#ifndef CALIBRATION_METHOD_H
#define CALIBRATION_METHOD_H

#include "predicter.h"

#include <complex>
#include <cstddef>
#include <vector>
#include <stdexcept>

class CalibrationMethod
{
	private:
		template<typename T, int NPol>
		class ValueArray
		{
			public:
				ValueArray(size_t nChannels, size_t nAntenna, size_t nTimesteps) :
					_nBaselines(nAntenna * (nAntenna-1) / 2),
					_nChannels(nChannels),
					_nAntenna(nAntenna),
					_nTimesteps(nTimesteps),
					_values(nChannels * _nBaselines * nTimesteps * NPol)
				{
				}
				
				T *ValuePtr(size_t antenna1, size_t antenna2, size_t timeIndex)
				{
					if(antenna1 == antenna2)
						throw std::runtime_error("Requested an auto-correlation");
					if(antenna1 >= _nAntenna || antenna2 >= _nAntenna || timeIndex >= _nTimesteps)
						throw std::runtime_error("Requested index value out of bounds");
					if(antenna1 > antenna2)
						std::swap(antenna1, antenna2);
					return &_values[NPol * _nChannels * (BaselineIndex(antenna1, antenna2) + timeIndex * _nBaselines)];
				}
				void SetAll(const T &newValue)
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
				std::vector<T> _values;
		};
		
		typedef ValueArray<std::complex<double>, 4> DataArray;
		typedef ValueArray<double, 1> WeightArray;
		
	public:
		CalibrationMethod(size_t nChannels, size_t nAntenna, size_t nTimesteps);
		
		void AddData(const std::complex<float> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2, size_t timestep);
		
		size_t Execute(double precisionLimit, size_t nIter);
		
		const std::complex<double> &JonesSolution(size_t antenna, size_t channel, size_t polarization) const
		{
			return _jonesSolutions[(antenna * _nChannels + channel) * 4 + polarization];
		}
		void SolutionSingularValue(size_t antenna, size_t channel, double &s1, double &s2) const
		{
			singularValues2x2(&_jonesSolutions[(antenna * _nChannels + channel) * 4], s1, s2);
		}

		static std::string MatrixToString(const std::complex<double> *matrix);
	private:
		void calculateNextIter(size_t ant, std::complex<double> *nextJones);
		
		void applyWeightsToData();
		double totalDistance(size_t antenna);
		double totalDistance();
		void reportDistances();
		
		static bool multiplyWithInverse2x2(std::complex<double>* lhs, const std::complex<double>* rhs);
		static void singularValues2x2(const std::complex<double>* matrix, double &e1, double &e2);
		
		DataArray _data, _model;
		WeightArray _weights, _weightSums;
		std::vector<std::complex<double> > _jonesSolutions;
		size_t _nChannels, _nAntenna, _nTimesteps;
		
		static void aTimesB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
		{
			dest[0] = lhs[0] * rhs[0] + lhs[1] * rhs[2];
			dest[1] = lhs[0] * rhs[1] + lhs[1] * rhs[3];
			dest[2] = lhs[2] * rhs[0] + lhs[3] * rhs[2];
			dest[3] = lhs[2] * rhs[1] + lhs[3] * rhs[3];
		}
		
		static void plusATimesB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
		{
			dest[0] += lhs[0] * rhs[0] + lhs[1] * rhs[2];
			dest[1] += lhs[0] * rhs[1] + lhs[1] * rhs[3];
			dest[2] += lhs[2] * rhs[0] + lhs[3] * rhs[2];
			dest[3] += lhs[2] * rhs[1] + lhs[3] * rhs[3];
		}
		
		static void aTimesHermB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
		{
			dest[0] = lhs[0] * std::conj(rhs[0]) + lhs[1] * std::conj(rhs[1]);
			dest[1] = lhs[0] * std::conj(rhs[2]) + lhs[1] * std::conj(rhs[3]);
			dest[2] = lhs[2] * std::conj(rhs[0]) + lhs[3] * std::conj(rhs[1]);
			dest[3] = lhs[2] * std::conj(rhs[2]) + lhs[3] * std::conj(rhs[3]);
		}

		static void hermATimesHermB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
		{
			dest[0] = std::conj(rhs[0]) * std::conj(lhs[0]) + std::conj(rhs[2]) * std::conj(lhs[1]);
			dest[1] = std::conj(rhs[0]) * std::conj(lhs[2]) + std::conj(rhs[2]) * std::conj(lhs[3]);
			dest[2] = std::conj(rhs[1]) * std::conj(lhs[0]) + std::conj(rhs[3]) * std::conj(lhs[1]);
			dest[3] = std::conj(rhs[1]) * std::conj(lhs[2]) + std::conj(rhs[3]) * std::conj(lhs[3]);
		}
		
		static void plusHermATimesB(std::complex<double> *dest, const std::complex<double> *lhs, const std::complex<double> *rhs)
		{
			dest[0] += std::conj(lhs[0]) * rhs[0] + std::conj(lhs[2]) * rhs[2];
			dest[1] += std::conj(lhs[0]) * rhs[1] + std::conj(lhs[2]) * rhs[3];
			dest[2] += std::conj(lhs[1]) * rhs[0] + std::conj(lhs[3]) * rhs[2];
			dest[3] += std::conj(lhs[1]) * rhs[1] + std::conj(lhs[3]) * rhs[3];
		}
};

#endif
