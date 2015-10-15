#include <iostream>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include "banddata.h"
#include "multibanddata.h"

int main(int argc, char* argv[])
{
	if(argc != 5)
	{
		std::cout << "Syntax: nwlayers <ms> <width> <height> <scale>\n";
		return -1;
	}
	casacore::MeasurementSet ms(argv[1]);
	const size_t imageWidth = atoi(argv[2]), imageHeight = atoi(argv[3]);
	double pixelSize = atof(argv[4])*(M_PI/180.0);
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casacore::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	
	const MultiBandData partBandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
	
	casacore::ROScalarColumn<int> ant1Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA1));
	casacore::ROScalarColumn<int> ant2Column(ms, ms.columnName(casacore::MSMainEnums::ANTENNA2));
	casacore::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casacore::MSMainEnums::DATA_DESC_ID));
	casacore::ROScalarColumn<int> fieldIdColumn(ms, ms.columnName(casacore::MSMainEnums::FIELD_ID));
	casacore::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casacore::MSMainEnums::FLAG));
	casacore::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casacore::MSMainEnums::UVW));
	
	double maxW= -1e100, minW = 1e100;
	double maxBaseline = 0.0;
	casacore::IPosition flagShape(flagColumn.shape(0));
	casacore::Array<bool> flagArray(flagShape);
	size_t polarizationCount = flagShape[0];
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			const BandData& curBand = partBandData[dataDescIdColumn(row)];
			casacore::Vector<double> uvwArray = uvwColumn(row);
			double uInM = uvwArray(0), vInM = uvwArray(1), wInM = uvwArray(2);
			double wHi = fabs(wInM / curBand.SmallestWavelength());
			double wLo = fabs(wInM / curBand.LongestWavelength());
			double baselineInM = sqrt(uInM*uInM + vInM*vInM + wInM*wInM);
			if(wHi > maxW || wLo < minW || baselineInM / curBand.SmallestWavelength() > maxBaseline)
			{
				flagColumn.get(row, flagArray);
				const bool* flagArrayData = flagArray.cbegin();
				for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
				{
					if(!*flagArrayData)
					{
						const double wavelength = curBand.ChannelWavelength(ch);
						maxW = std::max(maxW, fabs(wInM / wavelength));
						minW = std::min(minW, fabs(wInM / wavelength));
						maxBaseline = std::max(maxBaseline, baselineInM / wavelength);
					}
					flagArrayData += polarizationCount;
				}
			}
		}
	}

	double
		maxL = imageWidth * pixelSize * 0.5,
		maxM = imageHeight * pixelSize * 0.5,
		lmSq = maxL * maxL + maxM * maxM;
	double radiansForAllLayers;
	if(lmSq < 1.0)
		radiansForAllLayers = 2 * M_PI * (maxW - minW) * (1.0 - sqrt(1.0 - lmSq));
	else
		radiansForAllLayers = 2 * M_PI * (maxW - minW);
	size_t suggestedGridSize = size_t(ceil(radiansForAllLayers));
	std::cout << suggestedGridSize << '\n';
}
	