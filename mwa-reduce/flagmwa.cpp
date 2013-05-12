#include <iostream>
#include <stdexcept>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: flagmwa <ms> [<timestepcount> [<sbcount> <centre-channels> <side-channels>]]\n"
			"Flags the start timestep and corrupted centre channels of each sub-band as needed for\n"
			"MWA. timestepcount defaults to 5, sbcount defaults to 24, centre and side channels to 1.\n";
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		size_t timestepCount = 5, sbCount=24, centreChannelCount=1, sideChannelCount=1;
		if(argc>2)
		{
			timestepCount = atoi(argv[2]);
			if(argc>3)
			{
				sbCount = atoi(argv[3]);
				if(argc>4)
				{
					centreChannelCount = atoi(argv[4]);
					if(argc > 5)
						sideChannelCount = atoi(argv[5]);
				}
			}
		}
		
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
		
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ROScalarColumn<double> timeColumn(ms, ms.columnName(MSMainEnums::TIME));
		ArrayColumn<bool> flagsColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition flagsShape = flagsColumn.shape(0);
		unsigned polarizationCount = flagsShape[0];
		
		size_t channelsPerSubband = channelCount / sbCount;
		if(channelsPerSubband*sbCount != channelCount)
			throw std::runtime_error("channels per subband != channelCount * subband count");
		std::cout << "Flagging..." << std::flush;
		
		/**
		 * Flag
		 */
		size_t flaggedCount = 0, totalCount = 0, timestepIndex = (size_t) -1;
		double curTime = -1;
		Array<bool> flags(flagsShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			if(timeColumn(rowIndex) != curTime)
			{
				curTime = timeColumn(rowIndex);
				timestepIndex++;
			}
			
			// Cross correlation?
			if(ant1Column.get(rowIndex) != ant2Column.get(rowIndex))
			{
				flagsColumn.get(rowIndex, flags);
				Array<bool>::contiter flagPtr = flags.cbegin();
				if(timestepIndex < timestepCount)
				{
					for(size_t ch=0; ch!=channelCount;++ch)
					{
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							*flagPtr = true;
							flaggedCount++;
							++flagPtr;
							++totalCount;
						}
					}
				} else {
					for(size_t sb=0; sb!=sbCount; ++sb)
					{
						for(size_t ch=0; ch!=channelsPerSubband; ++ch)
						{
							bool doFlag = false;
							if(std::abs((long) (ch*2-channelsPerSubband))*2 < (long) centreChannelCount) doFlag=true;
							if(ch < sideChannelCount || ch >= channelsPerSubband-sideChannelCount)
								doFlag = true;
							for(size_t p=0; p!=polarizationCount; ++p)
							{
								if(doFlag) { *flagPtr = true; flaggedCount++; }
								++flagPtr;
								++totalCount;
							}
						}
					}
				}
				flagsColumn.put(rowIndex, flags);
			}
		}
		
		std::cout << " DONE (Set " << (round(flaggedCount*1000.0/totalCount)/10.0) << "% flagged)\n";
	}
	
	return 0;
}
