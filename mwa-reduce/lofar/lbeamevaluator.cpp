#include "lbeamevaluator.h"

#include "../banddata.h"

#include <measures/TableMeasures/ArrayMeasColumn.h>

#include <StationResponse/LofarMetaDataUtil.h>

using namespace LOFAR::StationResponse;

LBeamEvaluator::LBeamEvaluator(casa::MeasurementSet& ms) : _ms(ms)
{
	BandData band(_ms.spectralWindow());
	_subbandFrequency = band.CentreFrequency();
	
	casa::MSAntenna aTable(ms.antenna());
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	_arrayPos = antPosColumn(0);
	
	casa::MSField fieldTable(ms.field());
	casa::ROScalarMeasColumn<casa::MDirection> delayDirColumn(fieldTable, casa::MSField::columnName(casa::MSFieldEnums::DELAY_DIR));
	if(fieldTable.nrow() != 1)
		throw std::runtime_error("Set has multiple fields");
	_delayDir = delayDirColumn(0);
	
	if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
		casa::ROArrayMeasColumn<casa::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
		_tileBeamDir = *(tileBeamDirColumn(0).data());
	} else {
		throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
	}
	
	_j2000ToITRFRef = casa::MDirection::Convert(_j2000Ref, _itrfRef);
	
	_stations.resize(aTable.nrow());
	readStations(ms, _stations.begin());
}

LBeamEvaluator::~LBeamEvaluator()
{ }

void LBeamEvaluator::SetTime(const casa::MEpoch& time)
{
	_time = time;
	_timeAsDouble = _time.getValue().get()*86400.0;	
	
	_frame = casa::MeasFrame(_arrayPos, _time);
	_j2000Ref = casa::MDirection::Ref(casa::MDirection::J2000, _frame);
	_itrfRef = casa::MDirection::Ref(casa::MDirection::ITRF, _frame);
	_j2000ToITRFRef = casa::MDirection::Convert(_j2000Ref, _itrfRef);
	dirToITRF(_delayDir, _station0);
	dirToITRF(_tileBeamDir, _tile0);
}

void LBeamEvaluator::Evaluate(double ra, double dec, double frequency, size_t antennaIndex, MC2x2& beamValues)
{
	static const casa::Unit radUnit("rad");
	casa::MDirection imageDir(casa::MVDirection(
		casa::Quantity(ra, radUnit),
		casa::Quantity(dec,radUnit)),
		_j2000Ref);

	vector3r_t itrfDirection;
	dirToITRF(imageDir, itrfDirection);
	
	matrix22c_t gainMatrix = _stations[antennaIndex]->response(_timeAsDouble, frequency, itrfDirection, _subbandFrequency, _station0, _tile0);
	beamValues.Data()[0] = gainMatrix[0][0];
	beamValues.Data()[1] = gainMatrix[0][1];
	beamValues.Data()[2] = gainMatrix[1][0];
	beamValues.Data()[3] = gainMatrix[1][1];
}

void LBeamEvaluator::Evaluate(const LBeamEvaluator::PrecalcPosInfo& posInfo, double frequency, size_t antennaIndex, MC2x2& beamValues)
{
	matrix22c_t gainMatrix = _stations[antennaIndex]->response(_timeAsDouble, frequency, posInfo.itrfDirection, _subbandFrequency, _station0, _tile0);
	beamValues.Data()[0] = gainMatrix[0][0];
	beamValues.Data()[1] = gainMatrix[0][1];
	beamValues.Data()[2] = gainMatrix[1][0];
	beamValues.Data()[3] = gainMatrix[1][1];
}

void LBeamEvaluator::EvaluateFullArray(const LBeamEvaluator::PrecalcPosInfo& posInfo, double frequency, MC2x2& beamValues)
{
	beamValues = MC2x2::Zero();
	for(Station::Ptr s : _stations)
	{
		matrix22c_t gainMatrix = s->response(_timeAsDouble, frequency, posInfo.itrfDirection, _subbandFrequency, _station0, _tile0);
		beamValues.Data()[0] += gainMatrix[0][0];
		beamValues.Data()[1] += gainMatrix[0][1];
		beamValues.Data()[2] += gainMatrix[1][0];
		beamValues.Data()[3] += gainMatrix[1][1];
	}
	beamValues /= double(_stations.size());
}

void LBeamEvaluator::PrecalculatePositionInfo(LBeamEvaluator::PrecalcPosInfo& posInfo, double raRad, double decRad)
{
	static const casa::Unit radUnit("rad");
	casa::MDirection imageDir(casa::MVDirection(
		casa::Quantity(raRad, radUnit),
		casa::Quantity(decRad,radUnit)),
		_j2000Ref);

	dirToITRF(imageDir, posInfo.itrfDirection);
}
