#ifndef LBEAM_EVALUATOR_H
#define LBEAM_EVALUATOR_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include <measures/Measures/MEpoch.h>

#include "../matrix2x2.h"

#ifdef HAVE_LOFAR_BEAM
#include <StationResponse/Station.h>
#endif

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>

#include <memory>

class LBeamEvaluator
{
public:
	class PrecalcPosInfo
	{
	public:
		friend class LBeamEvaluator;
	private:
#ifdef HAVE_LOFAR_BEAM
		LOFAR::StationResponse::vector3r_t itrfDirection;
#endif
	};
	
	LBeamEvaluator(casa::MeasurementSet& ms);
	~LBeamEvaluator();

	void Evaluate(double ra, double dec, double frequency, size_t antennaIndex, MC2x2& beamValues);
		
	void Evaluate(const PrecalcPosInfo& posInfo, double frequency, size_t antennaIndex, MC2x2& beamValues);
	
	void EvaluateFullArray(const PrecalcPosInfo& posInfo, double frequency, MC2x2& beamValues);
		
	void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, double raRad, double decRad);
	
	void SetTime(const casa::MEpoch& time);
	
	const casa::MEpoch& Time() const { return _time; }
	
private:
	casa::MeasurementSet _ms;
	casa::MEpoch _time;
	double _timeAsDouble;
	
#ifdef HAVE_LOFAR_BEAM
	std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	double _subbandFrequency;
	casa::MDirection _delayDir, _tileBeamDir;
	casa::MPosition _arrayPos;
	casa::MeasFrame _frame;
	casa::MDirection::Ref _j2000Ref, _itrfRef;
	casa::MDirection::Convert _j2000ToITRFRef;
	LOFAR::StationResponse::vector3r_t _station0, _tile0;

	void dirToITRF(const casa::MDirection& dir, LOFAR::StationResponse::vector3r_t& itrf)
	{
		casa::MDirection itrfDir = _j2000ToITRFRef(dir);
		casa::Vector<double> itrfVal = itrfDir.getValue().getValue();
		itrf[0] = itrfVal[0];
		itrf[1] = itrfVal[1];
		itrf[2] = itrfVal[2];
	}
#endif

};

#ifndef HAVE_LOFAR_BEAM

// If the LOFAR beam library is not available, replace methods by dummies
inline LBeamEvaluator::LBeamEvaluator(casa::MeasurementSet&) { }

inline LBeamEvaluator::~LBeamEvaluator() { }

inline void LBeamEvaluator::Evaluate(double ra, double dec, double f, size_t a, MC2x2& beamValues)
{
	beamValues = MC2x2::Unity();
}
		
inline void LBeamEvaluator::Evaluate(const PrecalcPosInfo& p, double f, size_t a, MC2x2& beamValues)
{
	beamValues = MC2x2::Unity();
}
	
inline void LBeamEvaluator::PrecalculatePositionInfo(PrecalcPosInfo& p, double ra, double dec) { }
	
inline void LBeamEvaluator::SetTime(const casa::MEpoch& time) { _time = time; }

#endif // HAVE_LOFAR_BEAM

#endif
