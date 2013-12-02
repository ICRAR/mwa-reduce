#ifndef PREDICTER_H
#define PREDICTER_H

#include <complex>
#include <vector>

class Predicter
{
	public:
		typedef double NumType;
		typedef std::complex<NumType> CNumType;
		
		Predicter(NumType phaseCentreRA, NumType phaseCentreDec, NumType startFrequency, NumType endFrequency, size_t channelCount) :
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
		void Initialize(class ModelSource &source, class BeamEvaluator *beamEvaluator = 0);
		void Initialize(class Model &model, const std::string &solutionFile = std::string(), class BeamEvaluator *beamEvaluator = 0);
		void ReportSources(class Model& model);
		void UpdateBeam(class Model& model);
		
		//CNumType Predict(const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
		//CNumType Predict(const class Model &model, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex);
		
		void Predict4(CNumType *dest, const class ModelSource &source, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2);
		void Predict4(CNumType *dest, const class Model &model, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2);
		
		NumType TotalFlux(size_t p) { return std::fabs(_totalFlux[p]); }
	private:
		void initialize(class ModelComponent &source);
		void updateBeam(class ModelComponent &source);
		
		void predict4(CNumType *dest, const class ModelComponent& component, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2);
		struct SourceParameters
		{
			NumType l, m, lmsqrt, *brightness;
			CNumType *beamValues;
		};
		void readSolutions(const std::string& solutionFile);
		
		NumType _ra0, _dec0, _startFrequency, _endFrequency;
		size_t _channelCount;
		class BeamEvaluator *_beamEvaluator;
		CNumType _totalFlux[4];
		std::vector<std::complex<double>> _rhsSolutions;
};

#endif
