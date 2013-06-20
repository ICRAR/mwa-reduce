#include "model.h"
#include "tilebeam.h"
#include "banddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cout << "Usage: modelbeam <ms> <model-in> <model-out>";
		return -1;
	}
	
	/**
		* Read some meta data from the measurement set
		*/
	casa::MeasurementSet ms(argv[1]);
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	casa::MPosition ant1Pos = antPosColumn(0);
	
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	casa::MEpoch time = timeColumn(0);
	
	casa::MSField fTable(ms.field());
	if(fTable.nrow() != 1) throw std::runtime_error("Need exactly one field in set");
	casa::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	casa::MDirection refDir = refDirColumn(0);
	
	Model model(argv[2]);
	BandData band(ms.spectralWindow());
	double frequency = (band.BandStart()+ band.BandEnd()) * 0.5;
	
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
	TileBeam tileBeam(delays);
		
	for(Model::iterator sourceIter=model.begin(); sourceIter!=model.end(); ++sourceIter)
	{
		double ra = sourceIter->PosRA();
		double dec = sourceIter->PosDec();
		SpectralEnergyDistribution &sed = sourceIter->SED();
		
		double x, y;
		tileBeam.AnalyticGain(refDir, time, ant1Pos, ra, dec, frequency, x, y);
		sed.SetConstantBeam(x, 0.0, 0.0, y);
	}
	
	model.Save(argv[3]);
}
