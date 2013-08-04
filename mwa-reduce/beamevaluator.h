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
		void EvaluateAbsGain(double ra, double dec, double frequency, NumType* absGain)
		{
			std::complex<double> gains[4];
			EvaluateGain(ra, dec, frequency, gains);
			for(size_t p=0; p!=4; ++p)
				absGain[p] = std::fabs(gains[p]);
		}
		
		template<typename NumType>
		void AbsToApparent(double ra, double dec, double frequency, NumType* pixelValues)
		{
			double gains[4];
			EvaluateAbsGain(ra, dec, frequency, gains);
			
			// Calculate A D A^H
			
			NumType part[4];
			part[0] = gains[0] * pixelValues[0] + gains[1] * pixelValues[2];
			part[1] = gains[0] * pixelValues[1] + gains[1] * pixelValues[3];
			part[2] = gains[2] * pixelValues[0] + gains[3] * pixelValues[2];
			part[3] = gains[2] * pixelValues[1] + gains[3] * pixelValues[3];

			NumType temp[4];
			temp[0] = part[0] * gains[0] + part[1] * gains[1];
			temp[1] = part[0] * gains[2] + part[1] * gains[3];
			temp[2] = part[2] * gains[0] + part[3] * gains[1];
			temp[3] = part[2] * gains[2] + part[3] * gains[3];
			
			for(size_t p=0; p!=4; ++p)
				pixelValues[p] = temp[p];
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, double frequency, NumType* data)
		{
			double gains[4];
			EvaluateAbsGain(ra, dec, frequency, gains);
			
			// Calculate A^H^1 D A^1
			
			double overDeterminant = 1.0 / (gains[0]*gains[3] - gains[1]*gains[2]);
			double invGain[4];
			invGain[0] = overDeterminant * gains[3];
			invGain[1] = -overDeterminant * gains[1];
			invGain[2] = -overDeterminant * gains[2];
			invGain[3] = overDeterminant * gains[0];
			
			double part[4];
			part[0] = invGain[0] * data[0] + invGain[1] * data[2];
			part[1] = invGain[0] * data[1] + invGain[1] * data[3];
			part[2] = invGain[2] * data[0] + invGain[3] * data[2];
			part[3] = invGain[2] * data[1] + invGain[3] * data[3];

			double temp[4];
			temp[0] = part[0] * invGain[0] + part[1] * invGain[1];
			temp[1] = part[0] * invGain[2] + part[1] * invGain[3];
			temp[2] = part[2] * invGain[0] + part[3] * invGain[1];
			temp[3] = part[2] * invGain[2] + part[3] * invGain[3];
			
			for(size_t p=0; p!=4; ++p)
				data[p] = temp[p];
		}
	private:
		std::unique_ptr<TileBeam> _tileBeam;
		casa::MPosition _ant1Pos;
		casa::MDirection _refDir;
		casa::MEpoch _time;
		double _frequency;
};

#endif
