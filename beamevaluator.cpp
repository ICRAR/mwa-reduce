#include <stdexcept>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/tables/Tables/TableRecord.h>

#include "beamevaluator.h"
#include "banddata.h"

BeamEvaluator::BeamEvaluator(casacore::MeasurementSet& ms, bool reportDelays, const std::string& mwaPath)
{
	/**
		* Read some meta data from the measurement set
		*/
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casacore::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	
	casacore::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	_ant1Pos = antPosColumn(0);
	
	casacore::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	_time = timeColumn(0);
	
	BandData band(ms.spectralWindow());
	_frequency = (band.BandStart()+ band.BandEnd()) * 0.5;
	
	casacore::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
	casacore::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
	casacore::Array<int> delaysArr = delaysCol(0);
	casacore::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
	double delays[16];
	if(reportDelays)
		std::cout << "Delays: [";
	for(int i=0; i!=16; ++i)
	{
		delays[i] = delaysArrPtr[i];
		if(reportDelays)
		{
			std::cout << delays[i];
			if(i != 15) std::cout << ',';
		}
	}
	if(reportDelays)
		std::cout << "]\n";
	_tileBeam.reset(new TileBeam(delays, false, mwaPath));
}

void BeamEvaluator::EvaluateAbsToApparentGain(double ra, double dec, double frequency, std::complex<double>* gains)
{
	_tileBeam->ArrayResponse(_time, _ant1Pos, ra, dec, frequency, gains);
}

#ifdef CUDA_SUPPORT
void BeamEvaluator::EvaluateAbsToApparentGain(const std::vector<PrecalcPosInfo>& posInfo, std::complex<double> *gains, size_t startChannel, size_t endChannel, size_t channelCount, size_t startFrequency, size_t endFrequency){
	
	size_t n_components = posInfo.size();
	double *az = new double[n_components];
	double *ze = new double[n_components];
	for(size_t i = 0; i < n_components; i++){
		az[i] = posInfo[i].azimuth;
		ze[i] = posInfo[i].zenithAngle;
	}
	_tileBeam->ArrayResponse(ze, az, n_components, gains, startChannel, endChannel, channelCount, startFrequency, endFrequency);
	delete[] az;
	delete[] ze;
}
#endif