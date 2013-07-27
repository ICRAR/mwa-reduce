#ifndef BEAM_EVALUATOR_H
#define BEAM_EVALUATOR_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>

#include <memory>

#include "tilebeam.h"

class BeamEvaluator
{
	public:
		BeamEvaluator(casa::MeasurementSet &ms);
		
		void EvaluateGain(double ra, double dec, std::complex<double> *gains)
		{
			EvaluateGain(ra, dec, _frequency, gains);
		}
		
		void EvaluateGain(double ra, double dec, double *gains)
		{
			EvaluateGain(ra, dec, _frequency, gains);
		}
		
		void EvaluateGain(double ra, double dec, double frequency, std::complex<double> *gains);
		
		void EvaluateGain(double ra, double dec, double frequency, double *gains)
		{
			std::complex<double> compGains[4];
			EvaluateGain(ra, dec, frequency, compGains);
			gains[0] = compGains[0].real();
			gains[1] = compGains[1].real();
			gains[2] = compGains[2].real();
			gains[3] = compGains[3].real();
		}
	private:
		std::unique_ptr<TileBeam> _tileBeam;
		casa::MPosition _ant1Pos;
		casa::MDirection _refDir;
		casa::MEpoch _time;
		double _frequency;
};

#endif
