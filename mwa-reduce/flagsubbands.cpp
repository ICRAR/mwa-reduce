#include <iostream>
#include <stdexcept>
#include <set>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cout << "Usage: flagsubbands <ms> <subband count> <list of subbands>\n"
			"The list of sub-bands are zero indexed and should be separated by spaces.\n"; 
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		size_t subbandCount = atoi(argv[2]);
		std::set<size_t> subbands;
		for(int argi=3;argi!=argc;++argi) subbands.insert(atoi(argv[argi]));
		
		/**
		 * Read some meta data from the measurement set
		 */
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(channelCount%subbandCount != 0) throw std::runtime_error("Channel counts is not divisable by subband count");
		size_t channelsPerSubband = channelCount / subbandCount;
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition flagsShape = flagColumn.shape(0);
		unsigned polarizationCount = flagsShape[0];
		
		/**
		 * Apply flags
		 */
		std::cout << "Flagging " << subbands.size() << " subbands... " << std::flush;
		Array<bool> flags(flagsShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			flagColumn.get(rowIndex, flags);
			Array<bool>::iterator flagPtr = flags.begin();
			for(size_t sb=0; sb!=subbandCount; ++sb)
			{
				bool isSubbandFlagged = (subbands.find(sb) != subbands.end());
				for(size_t ch=0; ch!=channelsPerSubband; ++ch)
				{
					for(size_t p=0; p!=polarizationCount; ++p)
					{
						if(isSubbandFlagged) *flagPtr = true;
						++flagPtr;
					}
				}
			}
			flagColumn.put(rowIndex, flags);
		}
		
		std::cout << "DONE\n";
	}
	
	return 0;
}
