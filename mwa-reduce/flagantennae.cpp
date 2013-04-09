#include <iostream>
#include <stdexcept>
#include <set>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: flagantennae <ms> <list of antennae>\n"
			"Flags all correlations in which the specified antennae contribute.\n"
			"The list of antennae is separated with spaces, and antennae are given by index, starting\n"
			"at zero.\n";
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		std::set<size_t> antennae;
		for(int argi=2;argi!=argc;++argi) antennae.insert(atoi(argv[argi]));
		
		/**
		 * Read some meta data from the measurement set
		 */
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<bool> flagsColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		IPosition flagsShape = flagsColumn.shape(0);
		unsigned polarizationCount = flagsShape[0];
		
		std::cout << "Flagging... " << std::flush;
		
		/**
		 * Flag
		 */
		Array<bool> flags(flagsShape);
		Array<bool>::contiter flagPtr = flags.cbegin();
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
			if(antennae.find(ant1Column.get(rowIndex)) != antennae.end() || antennae.find(ant2Column.get(rowIndex)) != antennae.end())
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
