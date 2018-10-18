#include "lbeamevaluator.h"

#ifdef HAVE_LOFAR_BEAM

#include "../banddata.h"

#include "../wsclean/logger.h"
#include "../units/radeccoord.h"

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>

#include <StationResponse/LofarMetaDataUtil.h>

using namespace LOFAR::StationResponse;

LBeamEvaluator::LBeamEvaluator(casacore::MeasurementSet& ms) : _ms(ms)
{
	BandData band(_ms.spectralWindow());
	_subbandFrequency = band.CentreFrequency();
	
	casacore::MSAntenna aTable(ms.antenna());
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	_arrayPos = antPosColumn(0);
	
	casacore::MSField fieldTable(ms.field());
	casacore::ROScalarMeasColumn<casacore::MDirection> delayDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
	if(fieldTable.nrow() != 1)
		throw std::runtime_error("Set has multiple fields");
	_delayDir = delayDirColumn(0);
	
	if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
		casacore::ROArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
		_tileBeamDir = *(tileBeamDirColumn(0).data());
	} else {
		_tileBeamDir = _delayDir;
	}
	
	Logger::Debug << "Using delay direction: " << RaDecCoord::RaDecToString(_tileBeamDir.getValue().getVector()[0], _tileBeamDir.getValue().getVector()[1]) << '\n';
	
	_stations.resize(aTable.nrow());
	readStations(ms, _stations.begin());
}

LBeamEvaluator::~LBeamEvaluator()
{ }

void LBeamEvaluator::SetTime(const casacore::MEpoch& time)
{
	_time = time;
	_timeAsDouble = _time.getValue().get()*86400.0;
	
	_frame = casacore::MeasFrame(_arrayPos, _time);
	_j2000Ref = casacore::MDirection::Ref(casacore::MDirection::J2000, _frame);
	_itrfRef = casacore::MDirection::Ref(casacore::MDirection::ITRF, _frame);
	_j2000ToITRFRef = casacore::MDirection::Convert(_j2000Ref, _itrfRef);
	dirToITRF(_delayDir, _station0);
	dirToITRF(_tileBeamDir, _tile0);
}

void LBeamEvaluator::Evaluate(double ra, double dec, double frequency, size_t antennaIndex, MC2x2& beamValues)
{
	static const casacore::Unit radUnit("rad");
	casacore::MDirection imageDir(casacore::MVDirection(
		casacore::Quantity(ra, radUnit),
		casacore::Quantity(dec,radUnit)),
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
	beamValues[0] = gainMatrix[0][0];
	beamValues[1] = gainMatrix[0][1];
	beamValues[2] = gainMatrix[1][0];
	beamValues[3] = gainMatrix[1][1];
}

void LBeamEvaluator::EvaluateFullArray(const LBeamEvaluator::PrecalcPosInfo& posInfo, double frequency, MC2x2& beamValues)
{
	beamValues = MC2x2::Zero();
	for(Station::Ptr s : _stations)
	{
		matrix22c_t gainMatrix = s->response(_timeAsDouble, frequency, posInfo.itrfDirection, _subbandFrequency, _station0, _tile0);
		beamValues[0] = beamValues[0] + gainMatrix[0][0];
		beamValues[1] = beamValues[1] + gainMatrix[0][1];
		beamValues[2] = beamValues[2] + gainMatrix[1][0];
		beamValues[3] = beamValues[3] + gainMatrix[1][1];
	}
	beamValues /= double(_stations.size());
}

void LBeamEvaluator::PrecalculatePositionInfo(LBeamEvaluator::PrecalcPosInfo& posInfo, double raRad, double decRad)
{
	static const casacore::Unit radUnit("rad");
	casacore::MDirection imageDir(casacore::MVDirection(
		casacore::Quantity(raRad, radUnit),
		casacore::Quantity(decRad,radUnit)),
		_j2000Ref);

	dirToITRF(imageDir, posInfo.itrfDirection);
}

void LBeamEvaluator::EvaluateFullCorrection(casacore::MeasurementSet& ms, double ra, double dec, const class BandData& band, MC2x2* beamValues)
{
	LBeamEvaluator evaluator(ms);
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	size_t nrow = ms.nrow();
	std::vector<size_t> count(band.ChannelCount(), 0);
	for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
		beamValues[ch] = MC2x2::Zero();
	
	for(size_t i=0; i!=nrow; ++i)
	{
		casacore::MEpoch time = timeColumn(i);
		if(time.getValue().get() != evaluator.Time().getValue().get())
		{
			evaluator.SetTime(time);
			LBeamEvaluator::PrecalcPosInfo posInfo;
			evaluator.PrecalculatePositionInfo(posInfo, ra, dec);
			for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
			{
				MC2x2 timeStepValue;
				evaluator.EvaluateFullArray(posInfo, band.ChannelFrequency(ch), timeStepValue);
				beamValues[ch] += timeStepValue;
				++count[ch];
			}
		}
	}
	for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
		beamValues[ch] /= double(count[ch]);
}


#endif
