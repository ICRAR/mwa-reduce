#ifndef FLUX_ACCUMULATOR_H
#define FLUX_ACCUMULATOR_H

#include "beamevaluator.h"
#include "matrix2x2.h"

class FluxAccumulator
{
public:
	FluxAccumulator(double l, double m, double lambda);
	void UpdateBeam(const std::complex<double>* beamGains, double ionG, double ionDL, double ionDM);
	void Add(const std::complex<double>* vis, const double visWeight, double u, double v, double w);
	
	double L() const { return _l; }
	double M() const { return _m; }
private:
	void accumulateBeforeBeamChange()
	{
		// Add the residual flux still left in the observation after ionpeel:
		// Calculate Flux += w B* V B  (from: w (B* B) B^-1 V B*^-1 (B* B))
		// w = data weight, B = beam weight, V = vis
		std::complex<double> temp[4], temp2[4];
		Matrix2x2::HermATimesB(temp, _beamGains, _accFluxesBeforeBeamChange);
		Matrix2x2::PlusATimesB(_accFluxes, temp, _beamGains);
		
		// Now add the flux subtracted by ionpeel
		// Flux += w * g * (B* B) ModelFlux (B* B)
		Matrix2x2::ATimesHermB(temp, _beamGains, _beamGains);
		Matrix2x2::ATimesB(temp2, temp, _modelFlux);
		Matrix2x2::ScalarMultiply(temp2, _accVisWeightBeforeBeamChange * _ionG);
		Matrix2x2::PlusATimesB(_accFluxes, temp2, temp);
		
		// Calculate Weight += w B* B
		Matrix2x2::HermATimesB(temp, _beamGains, _beamGains);
		Matrix2x2::HermATimesB(temp2, temp, temp);
		Matrix2x2::MultiplyAdd(_accWeights, temp2, _accVisWeightBeforeBeamChange);
	}
	
	std::complex<double> _accFluxesBeforeBeamChange[4];
	std::complex<double> _accFluxes[4];
	std::complex<double> _accWeights[4];
	
	std::complex<double> _beamGains[4];
	std::complex<double> _modelFlux[4];
	
	double _accVisWeightBeforeBeamChange;
	
	double _ionG, _l, _m, _movedL, _movedM, _movedLMSqrt;
};

#endif
