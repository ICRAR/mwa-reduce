#ifndef BEAM_EVALUATOR_H
#define BEAM_EVALUATOR_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>

#include <memory>

#include "tilebeam.h"
#include "matrix2x2.h"

class BeamEvaluator
{
	public:
		BeamEvaluator(casa::MeasurementSet &ms);
		
		void EvaluateApparentToAbsGain(double ra, double dec, std::complex<double> *gains)
		{
			EvaluateApparentToAbsGain(ra, dec, _frequency, gains);
		}
		
		void EvaluateAbsToApparentGain(double ra, double dec, std::complex<double> *gains)
		{
			EvaluateAbsToApparentGain(ra, dec, _frequency, gains);
		}
		
		void EvaluateApparentToAbsGain(double ra, double dec, double frequency, std::complex<double> *gains)
		{
			EvaluateAbsToApparentGain(ra, dec, frequency, gains);
			Matrix2x2::Invert(gains);
		}
		
		void EvaluateAbsToApparentGain(double ra, double dec, double frequency, std::complex<double> *gains);
		
		template<typename NumType>
		void AbsToApparent(double ra, double dec, double frequency, std::complex<NumType>* pixelValues)
		{
			std::complex<NumType> gains[4], temp[4];
			EvaluateAbsToApparentGain(ra, dec, frequency, gains);
			
			// Calculate A D A^H
			Matrix2x2::ATimesB(temp, gains, pixelValues);
			Matrix2x2::ATimesHermB(pixelValues, temp, gains);
		}
		
		template<typename NumType>
		void AbsToApparent(double ra, double dec, double frequency, NumType* pixelValues)
		{
			std::complex<NumType> input[4], gains[4], temp[4];
			EvaluateAbsToApparentGain(ra, dec, frequency, gains);
			input[0] = pixelValues[0]; input[1] = pixelValues[1];
			input[2] = pixelValues[2]; input[3] = pixelValues[3];
			// Calculate A D A^H
			Matrix2x2::ATimesB(temp, gains, input);
			Matrix2x2::ATimesHermB(input, temp, gains);
			pixelValues[0] = input[0].real(); pixelValues[1] = input[1].real();
			pixelValues[2] = input[2].real(); pixelValues[3] = input[3].real();
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, double frequency, std::complex<NumType>* data)
		{
			std::complex<double> gains[4], temp[4];
			EvaluateApparentToAbsGain(ra, dec, frequency, gains);
			
			// Calculate A^1 D A^1^H
			Matrix2x2::ATimesB(temp, gains, data);
			Matrix2x2::ATimesHermB(data, temp, gains);
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, double frequency, NumType* pixelValues)
		{
			std::complex<double> input[4], gains[4], temp[4];
			EvaluateApparentToAbsGain(ra, dec, frequency, gains);
			input[0] = pixelValues[0]; input[1] = pixelValues[1];
			input[2] = pixelValues[2]; input[3] = pixelValues[3];
			// Calculate A D A^H
			Matrix2x2::ATimesB(temp, gains, input);
			Matrix2x2::ATimesHermB(input, temp, gains);
			pixelValues[0] = input[0].real(); pixelValues[1] = input[1].real();
			pixelValues[2] = input[2].real(); pixelValues[3] = input[3].real();
		}
		
	private:
		std::unique_ptr<TileBeam> _tileBeam;
		casa::MPosition _ant1Pos;
		casa::MEpoch _time;
		double _frequency;
};

#endif
