#ifndef LBEAM_EVALUATOR_H
#define LBEAM_EVALUATOR_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include <measures/Measures/MEpoch.h>

#include "../matrix2x2.h"

#include <StationResponse/Station.h>

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
		LOFAR::StationResponse::vector3r_t itrfDirection;
	};
	
	LBeamEvaluator(casa::MeasurementSet& ms);
	~LBeamEvaluator();

	void Evaluate(double ra, double dec, double frequency, size_t antennaIndex, MC2x2& beamValues);
		
	void Evaluate(const PrecalcPosInfo& posInfo, double frequency, size_t antennaIndex, MC2x2& beamValues);
	
	void PrecalculatePositionInfo(PrecalcPosInfo& posInfo, double raRad, double decRad);
	
	void SetTime(const casa::MEpoch& time);
	
private:
	casa::MeasurementSet _ms;
	casa::MEpoch _time;
	double _timeAsDouble;
	
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
};

#endif
