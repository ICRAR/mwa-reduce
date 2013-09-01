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
		
		template<typename DataFloatType>
		void AddData(const std::complex<DataFloatType> *data, const float* weights, const std::complex<double> *predictedValues, size_t antenna1, size_t antenna2, size_t timestep);
		
		void Execute(double& precisionLimit, size_t& nIter);
		
		const std::complex<double> &JonesSolution(size_t antenna, size_t channel, size_t polarization) const
		{
			return _jonesSolutions[(antenna * _nChannels + channel) * 4 + polarization];
		}
		
		void SolutionSingularValue(size_t antenna, size_t channel, double &s1, double &s2) const;

		static std::string MatrixToString(const std::complex<double> *matrix);
		
		void InitSolutionsToUnity();
		void InitSolutionsToNaN();
		void InitSolutions(const CalibrationMethod &source)
		{
			_jonesSolutions = source._jonesSolutions;
		}
		void SetOnlySolveDiag(bool onlySolveDiag) {
			_onlySolveDiag = onlySolveDiag;
		}
		void SetOnlySolveScalar(bool onlySolveScalar) {
			_onlySolveScalar = onlySolveScalar;
		}
		void SetOnlySolveRotation(bool onlySolveRotation) {
			_onlySolveRotation = onlySolveRotation;
		}
		bool OnlySolveRotation() const { return _onlySolveRotation; }
		
		static double DefaultMinAccuracy() { return 0.00001; }
		static double DefaultStoppingAccuracy() { return 0.0000001; }
	private:
		void calculateNextIter(size_t ant, std::complex<double> *nextJones);
		
		void applyWeightsToData();
		double totalDistance(size_t antenna);
		double totalDistance();
		void reportDistances();
		
		DataArray _data, _model;
		WeightArray _weights, _weightSums;
		std::vector<std::complex<double> > _jonesSolutions;
		size_t _nChannels, _nAntenna, _nTimesteps;
		bool _onlySolveDiag, _onlySolveScalar, _onlySolveRotation;
		
};

#endif
