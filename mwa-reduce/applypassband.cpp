#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: applypassband [-a] [-r] <passband-txt-file> <ms>\n-a will apply it on the auto-correlations instead of the default cross correlations.\n-r will reverse passband.\n";
	} else {
		bool onAutos = false, reverse = false;
		int argi = 1;
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-a")==0)
				onAutos = true;
			else if(strcmp(argv[argi], "-r")==0)
				reverse = true;
			else
				throw std::runtime_error("Could not parse a parameter");
			++argi;
		}
		std::ifstream passbandFile(argv[argi]); ++argi;
		MeasurementSet ms(argv[argi], Table::Update);
		
		/**
		 * Read some meta data from the measurement set
		 */
		//MSAntenna aTable = ms.antenna();
		//size_t antennaCount = aTable.nrow();
		
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		
		/**
		 * Read the pass band file
		 */
		std::vector<long double> correctionFacts[polarizationCount];
		while(passbandFile.good())
		{
			size_t rowIndex;
			passbandFile >> rowIndex;
			if(!passbandFile.good()) break;
			for(size_t p=0; p!=polarizationCount; ++p) {
				long double stddev;
				passbandFile >> stddev;
				if(onAutos) stddev = sqrt(stddev);
				if(reverse)
					correctionFacts[p].push_back(stddev);
				else
					correctionFacts[p].push_back(1.0/stddev);
			}
		}
		size_t
			channelsPerSubband = correctionFacts[0].size(),
			subBandCount = channelCount / channelsPerSubband;
		std::cout << "Passband file contains " << channelsPerSubband << " channels.\n";
		std::cout << "Measurement set contains " << channelCount << " channels.\n";
		std::cout << "Assuming set contains " << subBandCount << " subbands.\n";
		std::cout << "Correcting passband..." << std::flush;
		
		
		/**
		 * Apply corrections
		 */
		Array<complex_t> data(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			if((ant1Column.get(rowIndex) != ant2Column.get(rowIndex) && !onAutos) ||
				(ant1Column.get(rowIndex) == ant2Column.get(rowIndex) && onAutos))
			{
				dataColumn.get(rowIndex, data);
				Array<complex_t>::iterator dataPtr = data.begin();
				for(size_t sb=0; sb!=subBandCount; ++sb)
				{
					for(size_t ch=0; ch!=channelsPerSubband; ++ch)
					{
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							*dataPtr = complex_t(
								dataPtr->real() * correctionFacts[p][ch],
								dataPtr->imag() * correctionFacts[p][ch]
							);
							++dataPtr;
						}
					}
				}
				dataColumn.put(rowIndex, data);
			}
		}
		
		std::cout << " DONE\n";
	}
	
	return 0;
}
