#ifndef BEAM_EVALUATOR_H
#define BEAM_EVALUATOR_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>

#include <memory>

#include "tilebeam.h"
#include "matrix2x2.h"
#include "modelsource.h"
#include "model.h"

class BeamEvaluator
{
	public:
		typedef TileBeam::PrecalcPosInfo PrecalcPosInfo;
		
		BeamEvaluator(casa::MeasurementSet& ms, bool reportDelays = true);
		
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
		void AbsToApparent(double ra, double dec, std::complex<NumType>* pixelValues)
		{
			AbsToApparent<NumType>(ra, dec, _frequency, pixelValues);
		}
		
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
		void ApparentToAbs(double ra, double dec, std::complex<NumType>* data)
		{
			ApparentToAbs<NumType>(ra, dec, _frequency, data);
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, double frequency, std::complex<NumType>* data)
		{
			std::complex<double> gains[4], temp[4];
			EvaluateApparentToAbsGain(ra, dec, frequency, gains);
			
			Matrix2x2::ATimesB(temp, gains, data);
			Matrix2x2::ATimesHermB(data, temp, gains);
		}
		
		template<typename NumType>
		void ApparentToAbs(double ra, double dec, NumType* pixelValues)
		{
			ApparentToAbs<NumType>(ra, dec, _frequency, pixelValues);
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
		
		void SetTime(const casa::MEpoch& time) { _time = time; }
		const casa::MEpoch& Time() { return _time; }
		
		void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, double raRad, double decRad)
		{
			_tileBeam->PrecalculatePositionInfo(posInfo, _time, _ant1Pos, raRad, decRad);
		}
		
		void EvaluateApparentToAbsGain(const PrecalcPosInfo& posInfo, double frequency, std::complex<double> *gains)
		{
			EvaluateAbsToApparentGain(posInfo, frequency, gains);
			Matrix2x2::Invert(gains);
		}
		
		void EvaluateAbsToApparentGain(const PrecalcPosInfo& posInfo, double frequency, std::complex<double> *gains)
		{
			_tileBeam->AnalyticJones(posInfo, frequency, gains);
		}
	
		/*void AbsToApparent(ModelSource& source, double frequency)
		{
			for(ModelSource::iterator i=source.begin(); i!=source.end(); ++i)
			{
				ModelComponent& comp = *i;
				double vals[4];
				for(size_t p=0; p!=4; ++p)
					vals[p] = comp.SED().FluxAtFrequency(frequency, p);
				AbsToApparent(comp.PosRA(), comp.PosDec(), frequency, vals);
				comp.SetSED(SpectralEnergyDistribution(vals, frequency));
			}
		}
		
		void AbsToApparent(ModelSource& source)
		{
			AbsToApparent(source, _frequency);
		}
		
		void AbsToApparent(Model& model, double frequency)
		{
			for(Model::iterator i=model.begin(); i!=model.end(); ++i)
				AbsToApparent(*i, frequency);
		}
		
		void AbsToApparent(Model& model)
		{
			AbsToApparent(model, _frequency);
		}*/
	private:
		std::unique_ptr<TileBeam> _tileBeam;
		casa::MPosition _ant1Pos;
		casa::MEpoch _time;
		double _frequency;
};

#endif
