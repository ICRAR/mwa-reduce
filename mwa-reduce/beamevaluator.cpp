#include <stdexcept>

#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>

#include "beamevaluator.h"
#include "banddata.h"

BeamEvaluator::BeamEvaluator(casa::MeasurementSet& ms)
{
	/**
		* Read some meta data from the measurement set
		*/
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	_ant1Pos = antPosColumn(0);
	
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	_time = timeColumn(0);
	
	casa::MSField fTable(ms.field());
	if(fTable.nrow() != 1) throw std::runtime_error("Need exactly one field in set");
	casa::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	_refDir = refDirColumn(0);
	
	BandData band(ms.spectralWindow());
	_frequency = (band.BandStart()+ band.BandEnd()) * 0.5;
	
	casa::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
	casa::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
	casa::Array<int> delaysArr = delaysCol(0);
	casa::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
	double delays[16];
	std::cout << "Delays: [";
	for(int i=0; i!=16; ++i)
	{
		delays[i] = delaysArrPtr[i];
		std::cout << delays[i];
		if(i != 15) std::cout << ',';
	}
	std::cout << "]\n";
	_tileBeam.reset(new TileBeam(delays));
}

void BeamEvaluator::EvaluateGain(double ra, double dec, double frequency, std::complex<double>* gains)
{
	_tileBeam->AnalyticJones(_refDir, _time, _ant1Pos, ra, dec, frequency, gains);
}
