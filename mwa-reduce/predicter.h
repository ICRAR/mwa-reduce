#ifndef PREDICTER_H
#define PREDICTER_H

#include <complex>
#include <vector>

class Predicter
{
	public:
		typedef double NumType;
		typedef std::complex<NumType> CNumType;
		
		Predicter(NumType phaseCentreRA, NumType phaseCentreDec, NumType startFrequency, NumType endFrequency, size_t channelCount, bool applyBeam = false) :
			_ra0(phaseCentreRA), _dec0(phaseCentreDec), _startFrequency(startFrequency), _endFrequency(endFrequency), _channelCount(channelCount)
		{
			_totalFlux[0] = 0.0;
			_totalFlux[1] = 0.0;
			_totalFlux[2] = 0.0;
			_totalFlux[3] = 0.0;
		}
		
		/**
		 * Initializes the l and m position(s) of the source.
		 */
		void Initialize(class ModelSource &model, class BeamEvaluator *beamEvaluator = 0);
		void Initialize(class Model &model, const std::string &solutionFile = std::string(), class BeamEvaluator *beamEvaluator = 0);
		void ReportSources(class Model& model);
		
		CNumType Predict(const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
		CNumType Predict(const class Model &model, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
		
		void Predict4(CNumType *dest, const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2);
		void Predict4(CNumType *dest, const class Model &model, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2);
		
		double TotalFlux(size_t p) { return _totalFlux[p]; }
	private:
		void predict4(CNumType *dest, const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex);
		struct SourceParameters
		{
			NumType l, m, lmsqrt, *brightness;
			CNumType *beamValues;
		};
		void applyGain(double *dataVal, const std::complex<double> *gain);
		void readSolutions(const std::string& solutionFile);
		
		NumType _ra0, _dec0, _startFrequency, _endFrequency;
		size_t _channelCount;
		class BeamEvaluator *_beamEvaluator;
		double _totalFlux[4];
		std::vector<std::complex<double>> _solutions;
};

#endif
