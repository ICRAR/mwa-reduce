#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "banddata.h"
#include "sourcestrength.h"

using namespace casa;

int main(int argc, char **argv)
{
	if(argc != 5)
	{
		std::cout << "Usage: subtract <ms> <source-flux-density> <spectral index> <ref frequency>\n"
			"Subtracts the flux density from the real part of the visibilities. This 'peels' the\n"
			"source in the phase centre. Only affects cross-correlations.\n"
			"Frequency given in MHz, flux density in Jy.\n";
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		const long double
			sourceFluxDensity = atof(argv[2]),
			spectralIndex = atof(argv[3]),
			refFrequency = atof(argv[4]);
		const SourceStrength<long double>
			sourceStrength(sourceFluxDensity, spectralIndex, refFrequency*1000000.0);
		
		/**
		 * Read some meta data from the measurement set
		 */
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		/**
		 * Apply spectral index to source flux
		 */
		double sourceFlux[channelCount];
		for(size_t ch=0;ch!=channelCount;++ch)
			sourceFlux[ch] = sourceStrength.FluxAtFrequency(bandData.ChannelFrequency(ch));
		std::cout << "Subtracting source of " << sourceFlux[0] << " - " << sourceFlux[channelCount-1] << " Jy... " << std::flush;
		
		/**
		 * Subtract
		 */
		Array<complex_t> data(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				dataColumn.get(rowIndex, data);
				Array<complex_t>::iterator dataPtr = data.begin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					for(size_t p=0;p!=polarizationCount;++p)
					{
						if(polarizationCount!=4 || p==0 || p==3)
							*dataPtr = complex_t(dataPtr->real()-sourceFlux[ch], dataPtr->imag());
						++dataPtr;
					}
				}
				dataColumn.put(rowIndex, data);
			}
		}
		
		std::cout << "DONE\n";		
	}
}
