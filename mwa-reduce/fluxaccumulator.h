#ifndef FLUX_ACCUMULATOR_H
#define FLUX_ACCUMULATOR_H

#include "beamevaluator.h"
#include "matrix2x2.h"
#include "serializable.h"

class FluxAccumulator : public Serializable
{
public:
	/**
	 * This constructor is used when calculating flux from visibilities
	 */
	FluxAccumulator(double l, double m, double lambda, const std::complex<double>* modelFlux);
	
	/**
	 * This constructor is used when accumulating serialized files.
	 */
	FluxAccumulator(double lambda);
	
	void UpdateBeam(const std::complex<double>* beamGains, double ionG, double ionDL, double ionDM);
	void Add(const std::complex<double>* vis, const double visWeight, double u, double v, double w);
	
	double L() const { return _l; }
	double M() const { return _m; }
	
	void Finish()
	{
		accumulateBeforeBeamChange();
	}
	
	virtual void Serialize(std::ostream& stream) const
	{
		if(_accVisWeightBeforeBeamChange != 0.0)
			throw std::runtime_error("Serializing flux accumulator that has not been finished!");
		
		serializeMatrix(stream, _accFluxes);
		serializeMatrix(stream, _accWeights);
		
		serializeMatrix(stream, _modelFlux);
		
		SerializeToDouble(stream, _l);
		SerializeToDouble(stream, _m);
	}
	
	virtual void Unserialize(std::istream& stream)
	{
		_accVisWeightBeforeBeamChange = 0.0;
		for(size_t p=0; p!=4; ++p)
			_accFluxesBeforeBeamChange[p] = 0.0;
		_beamGains[0] = 1.0; _beamGains[1] = 0.0;
		_beamGains[2] = 0.0; _beamGains[3] = 1.0;
		
		unserializeMatrix(stream, _accFluxes);
		unserializeMatrix(stream, _accWeights);
		
		unserializeMatrix(stream, _modelFlux);
		
		_l = UnserializeDouble(stream);
		_m = UnserializeDouble(stream);
	}
	
	virtual void AccumulateFromStream(std::istream& stream)
	{
		unserializeAndAddMatrix(stream, _accFluxes);
		unserializeAndAddMatrix(stream, _accWeights);
		
		std::complex<double> model[4];
		unserializeMatrix(stream, model); // model
		
		UnserializeDouble(stream); // l
		UnserializeDouble(stream); // m
	}
	
	void GetFlux(double* stokesMatrix) const
	{
		// Calculate: W*FW sum(W*W)^-1
		// sum(W*W) is in variable _accWeights, and W*FW in _accFluxes.
		std::complex<double> weights[4];
		memcpy(weights, _accWeights, sizeof(std::complex<double>)*4);
		if(Matrix2x2::Invert(weights))
		{
			std::complex<double> correctedLinear[4];
			Matrix2x2::ATimesB(correctedLinear, weights, _accFluxes);
			Polarization::LinearToStokes(correctedLinear, stokesMatrix);
		}
		else {
			for(size_t p=0; p!=4; ++p)
				stokesMatrix[p] = std::numeric_limits<double>::quiet_NaN();
		}
	}
private:
	static void serializeMatrix(std::ostream& stream, const std::complex<double>* m)
	{
		for(size_t p=0; p!=4; ++p)
			SerializeToDoubleC(stream, m[p]);
	}
	static void unserializeMatrix(std::istream& stream, std::complex<double>* m)
	{
		for(size_t p=0; p!=4; ++p)
			m[p] = UnserializeDoubleC(stream);
	}
	static void unserializeAndAddMatrix(std::istream& stream, std::complex<double>* m)
	{
		for(size_t p=0; p!=4; ++p)
			m[p] += UnserializeDoubleC(stream);
	}
	
	void accumulateBeforeBeamChange()
	{
		if(std::isfinite(_ionG))
		{
			// Add the residual flux still left in the observation after ionpeel:
			// Calculate Flux += w B* V B  (from: w (B* B) B^-1 V B*^-1 (B* B))
			// w = data weight * g^2, B = beam weight, V = vis
			std::complex<double> temp[4], temp2[4];
			Matrix2x2::HermATimesB(temp, _beamGains, _accFluxesBeforeBeamChange);
			Matrix2x2::ScalarMultiply(temp, _ionG*_ionG);
			Matrix2x2::PlusATimesB(_accFluxes, temp, _beamGains);
			
			// Now add the flux subtracted by ionpeel
			// Flux += w * g * (B* B) ModelFlux (B* B)
			Matrix2x2::HermATimesB(temp, _beamGains, _beamGains);
			Matrix2x2::ATimesB(temp2, temp, _modelFlux);
			Matrix2x2::ScalarMultiply(temp2, _accVisWeightBeforeBeamChange * _ionG * (_ionG*_ionG));
			Matrix2x2::PlusATimesB(_accFluxes, temp2, temp);
			
			// Calculate Weight += w (B* B) (B* B)
			Matrix2x2::HermATimesB(temp, _beamGains, _beamGains);
			Matrix2x2::ATimesB(temp2, temp, temp); //Herm?!
			Matrix2x2::MultiplyAdd(_accWeights, temp2, _accVisWeightBeforeBeamChange * (_ionG*_ionG));
		}
		
		for(size_t p=0; p!=4; ++p)
			_accFluxesBeforeBeamChange[p] = std::complex<double>(0.0, 0.0);
		_accVisWeightBeforeBeamChange = 0.0;
	}
	
	std::complex<double> _accFluxesBeforeBeamChange[4];
	std::complex<double> _accFluxes[4];
	std::complex<double> _accWeights[4];
	
	std::complex<double> _beamGains[4];
	std::complex<double> _modelFlux[4];
	
	double _accVisWeightBeforeBeamChange;
	
	double _l, _m, _ionG, _movedL, _movedM, _movedLMSqrt;
};

#endif
