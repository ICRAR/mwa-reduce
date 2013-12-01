#ifndef CALIBRATION_METHOD_H
#define CALIBRATION_METHOD_H

#include "predicter.h"
#include "visibilityarray.h"

#include <complex>
#include <cstddef>
#include <vector>
#include <stdexcept>

class CalibrationMethod
{
	private:
		
		typedef VisibilityArray<std::complex<double>, 4> DataArray;
		typedef VisibilityArray<double, 1> WeightArray;
		
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
		
		static double DefaultMinAccuracy() { return 0.0001; }
		static double DefaultStoppingAccuracy() { return 0.000001; }
		static size_t DefaultNIter() { return 1000; }
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
