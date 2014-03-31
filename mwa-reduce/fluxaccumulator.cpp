#include "fluxaccumulator.h"

FluxAccumulator::FluxAccumulator(double l, double m, double lambda) :
	_accFluxesBeforeBeamChange({0.0, 0.0, 0.0, 0.0}),
	_accFluxes({0.0, 0.0, 0.0, 0.0}),
	_accWeights({0.0, 0.0, 0.0, 0.0}),
	_beamGains({0.0, 0.0, 0.0, 0.0}),
	_accVisWeightBeforeBeamChange(0.0)
{
}

void FluxAccumulator::UpdateBeam(const std::complex<double>* beamGains, double ionG, double ionDL, double ionDM)
{
	accumulateBeforeBeamChange();
	
	_movedL = _l+ionDL;
	_movedM = _m+ionDM;
	_movedLMSqrt = sqrt(1.0 - _movedL*_movedL - _movedM*_movedM)-1.0;
	_ionG = ionG;
	Matrix2x2::Assign(_beamGains, beamGains);
	for(size_t p=0; p!=4; ++p)
		_accFluxesBeforeBeamChange[p] = std::complex<double>(0.0, 0.0);
	_accVisWeightBeforeBeamChange = 0.0;
}

void FluxAccumulator::Add(const std::complex<double>* vis, const double visWeight, double u, double v, double w)
{
	double angle = 2.0*M_PI*(u*_movedL + v*_movedM + w*_movedLMSqrt);
	
	double sinAngle, cosAngle;
	sincos(angle, &sinAngle, &cosAngle);
	sinAngle *= visWeight;
	cosAngle *= visWeight;
	std::complex<double>
		fact(cosAngle - sinAngle, sinAngle + cosAngle),
		predict[4] = { fact, std::complex<double>(0.0), std::complex<double>(0.0), fact };
	Matrix2x2::PlusATimesHermB(_accFluxesBeforeBeamChange, vis, predict);
	Matrix2x2::PlusHermATimesB(_accFluxesBeforeBeamChange, vis, predict);
	// Weight is multiplied by factor of 2 because we add both normal and conjugate at once
	_accVisWeightBeforeBeamChange += 2.0 * visWeight;
}
