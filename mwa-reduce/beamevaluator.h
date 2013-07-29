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
		
		void EvaluateGain(double ra, double dec, double frequency, std::complex<double> *gains);
		
		template<typename NumType>
		void AbsToApparent(double ra, double dec, double frequency, NumType* pixelValues)
		{
			// TODO this should probably be done as "|G| d |G^H|" instead of d |GG&^H|
			// Doesn't matter for unpolarized sources (d is diagonal and commutative)
			std::complex<double> gains[4];
			EvaluateGain(ra, dec, frequency, gains);
			
			double gainSq[4];
			gainSq[0] = std::fabs(gains[0] * gains[0] + gains[1] * gains[1]);
			gainSq[1] = std::fabs(gains[0] * gains[2] + gains[1] * gains[3]);
			gainSq[2] = std::fabs(gains[2] * gains[0] + gains[3] * gains[1]);
			gainSq[3] = std::fabs(gains[2] * gains[2] + gains[3] * gains[3]);

			double temp[4];
			temp[0] = pixelValues[0] * gainSq[0] + pixelValues[1] * gainSq[2];
			temp[1] = pixelValues[0] * gainSq[1] + pixelValues[1] * gainSq[3];
			temp[2] = pixelValues[2] * gainSq[0] + pixelValues[3] * gainSq[2];
			temp[3] = pixelValues[2] * gainSq[1] + pixelValues[3] * gainSq[3];
			
			for(size_t p=0; p!=4; ++p)
				pixelValues[p] = temp[p];
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, double frequency, NumType* pixelValues)
		{
			std::complex<double> gains[4];
			EvaluateGain(ra, dec, frequency, gains);
			
			std::complex<double> overDeterminant = 1.0 / (gains[0]*gains[3] - gains[1]*gains[2]);
			std::complex<double> invGain[4];
			invGain[0] = overDeterminant * gains[3];
			invGain[1] = -overDeterminant * gains[1];
			invGain[2] = -overDeterminant * gains[2];
			invGain[3] = overDeterminant * gains[0];
			
			double gainSq[4];
			gainSq[0] = std::fabs(invGain[0] * invGain[0] + invGain[1] * invGain[1]);
			gainSq[1] = std::fabs(invGain[0] * invGain[2] + invGain[1] * invGain[3]);
			gainSq[2] = std::fabs(invGain[2] * invGain[0] + invGain[3] * invGain[1]);
			gainSq[3] = std::fabs(invGain[2] * invGain[2] + invGain[3] * invGain[3]);

			double temp[4];
			temp[0] = pixelValues[0] * gainSq[0] + pixelValues[1] * gainSq[2];
			temp[1] = pixelValues[0] * gainSq[1] + pixelValues[1] * gainSq[3];
			temp[2] = pixelValues[2] * gainSq[0] + pixelValues[3] * gainSq[2];
			temp[3] = pixelValues[2] * gainSq[1] + pixelValues[3] * gainSq[3];
			
			for(size_t p=0; p!=4; ++p)
				pixelValues[p] = temp[p];
		}
	private:
		std::unique_ptr<TileBeam> _tileBeam;
		casa::MPosition _ant1Pos;
		casa::MDirection _refDir;
		casa::MEpoch _time;
		double _frequency;
};

#endif
