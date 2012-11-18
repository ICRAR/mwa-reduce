#include <iostream>
#include <stdexcept>
#include <set>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "banddata.h"

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3 || (argc-2)%2 != 0)
	{
		std::cout << "Usage: flagantennae <ms> <list of baselines>\n"
			"Flags all given baselines.\n"
			"The list of baselines is separated with spaces, and antennae are given by index, starting\n"
			"at zero, e.g.: 'flagbaseline myms.ms 12 15 12 16' will flag the baselines 12x15 and 12x16.\n";
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		std::set<std::pair<size_t, size_t > > baselines;
		for(int argi=2;argi!=argc;argi+=2)
		{
			size_t a1 = atoi(argv[argi]), a2 = atoi(argv[argi+1]);
			if(a1 > a2) std::swap(a1, a2);
			baselines.insert(std::pair<size_t, size_t>(a1, a2));
		}
		
		/**
		 * Read some meta data from the measurement set
		 */
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<bool> flagsColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		IPosition flagsShape = flagsColumn.shape(0);
		unsigned polarizationCount = flagsShape[0];
		
		std::cout << "Flagging " << baselines.size() << " baselines... " << std::flush;
		
		/**
		 * Flag
		 */
		Array<bool> flags(flagsShape);
		Array<bool>::iterator flagPtr = flags.begin();
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			for(size_t p=0; p!=polarizationCount; ++p)
			{
				*flagPtr = true;
				++flagPtr;
			}
		}
		
		size_t crossCount = 0, autoCount = 0;
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Selected?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 > a2) std::swap(a1, a2);
			if(baselines.find(std::pair<size_t, size_t>(a1, a2)) != baselines.end())
			{
				if(ant1Column.get(rowIndex) == ant2Column.get(rowIndex))
					++autoCount;
				else
					++crossCount;
				flagsColumn.put(rowIndex, flags);
			}
		}
		
		std::cout << "DONE (selected " << crossCount << " cross- and " << autoCount << " auto-correlated timesteps)\n";
	}
	
	return 0;
}
