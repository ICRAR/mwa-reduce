#ifndef MULTIBANDDATA_H
#define MULTIBANDDATA_H

#include <stdexcept>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include "banddata.h"

/**
 * Contains information about a set of bands. This follows the CASA Measurement
 * Set model; one MultiBandData instance can contain the band information contained
 * in the CASA Measurement Set.
 */
class MultiBandData
{
	public:
		using iterator = std::vector<BandData>::iterator;
		using const_iterator = std::vector<BandData>::const_iterator;
		
		/**
		 * Construct an empty MultiBandData.
		 */
		MultiBandData()
		{ }
		
		/**
		 * Construct a MultiBandData from the Measurement Set tables.
		 * @param spwTable The spectral window table of a measurement set.
		 * @param dataDescTable The data description table of a measurement set.
		 */
		MultiBandData(casacore::MSSpectralWindow& spwTable, casacore::MSDataDescription& dataDescTable);
		
		/**
		 * Construct a MultiBandData from another instance but only select a part of each
		 * band data.
		 * @param source Other instance that will be partially copied.
		 * @param startChannel Start of channel range to initialize this instance with.
		 * @param endChannel End of channel range (exclusive) to initialize this instance with.
		 */
		MultiBandData(const MultiBandData& source, size_t startChannel, size_t endChannel);
		
		/**
		 * Index operator to retrieve a band data given a dataDescID.
		 * @param dataDescID A valid data description ID for which the band is returned.
		 * @returns The BandData for the requested band.
		 */
		const BandData& operator[](size_t dataDescID) const
		{
			return _bandData[_dataDescToBand[dataDescID]];
		}
		
		/**
		 * Retrieve the first band.
		 * @returns The first band.
		 */
		const BandData& FirstBand() const
		{
			return _bandData.front();
		}
		
		/**
		 * Get number of bands stored.
		 * @returns Number of bands.
		 */
		size_t BandCount() const
		{
			return _bandData.size();
		}
		
		/**
		 * Returns the unique number of data description IDs.
		 * @returns Unique number of data desc IDs.
		 */
		size_t DataDescCount() const
		{
			return _dataDescToBand.size();
		}
		
		/**
		 * Get lowest frequency.
		 * @returns The channel frequency of the channel with lowest frequency.
		 */
		double LowestFrequency() const
		{
			double freq = _bandData[0].LowestFrequency();
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::min(freq, _bandData[i].LowestFrequency());
			return freq;
		}
		
		/**
		 * Get centre frequency.
		 * @returns (LowestFrequency() + HighestFrequency()) * 0.5.
		 */
		double CentreFrequency() const
		{
			return (LowestFrequency() + HighestFrequency()) * 0.5;
		}
		
		/**
		 * Get highest frequency.
		 * @returns The channel frequency of the channel with highest frequency.
		 */
		double HighestFrequency() const
		{
			double freq = _bandData[0].HighestFrequency();
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::max(freq, _bandData[i].HighestFrequency());
			return freq;
		}
		
		/**
		 * Get total bandwidth covered.
		 * @returns BandEnd() - BandStart().
		 */
		double Bandwidth() const
		{
			return BandEnd() - BandStart();
		}
		
		/**
		 * Get the start frequency of the lowest frequency channel.
		 * @return Start of covered bandwidth.
		 */
		double BandStart() const
		{
			double freq = std::min(_bandData[0].BandStart(), _bandData[0].BandEnd());
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::min(freq, std::min(_bandData[i].BandStart(), _bandData[i].BandEnd()));
			return freq;
		}
		
		/**
		 * Get the end frequency of the highest frequency channel.
		 * @return End of covered bandwidth.
		 */
		double BandEnd() const
		{
			double freq = std::max(_bandData[0].BandStart(), _bandData[0].BandEnd());
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::max(freq, std::max(_bandData[i].BandStart(), _bandData[i].BandEnd()));
			return freq;
		}
		
		/**
		 * Get the maximum number of channels in a band.
		 * @returns Maximum number of channels.
		 */
		size_t MaxChannels() const
		{
			size_t maxChannels = 0;
			for(std::vector<BandData>::const_iterator i=_bandData.begin(); i!=_bandData.end(); ++i)
			{
				if(i->ChannelCount() > maxChannels)
					maxChannels = i->ChannelCount();
			}
			return maxChannels;
		}
		
		/**
		 * Map a dataDescId to the corresponding band index.
		 * @param dataDescId A dataDescId as e.g. used in a main table.
		 * @returns The band index, which is equal to the row index in the spw
		 * table that describes the band in a measurement set.
		 */
		size_t GetBandIndex(size_t dataDescId) const { return _dataDescToBand[dataDescId]; }
		
		/**
		 * Compose a list of dataDescIds that are used in the measurement set.
		 * "Used" here means it is references from the main table.
		 * @param mainTable the measurement set.
		 * @returns Set of used dataDescIds.
		 */
		std::set<size_t> GetUsedDataDescIds(casacore::MeasurementSet& mainTable) const;
		
		iterator begin() { return _bandData.begin(); }
		const_iterator begin() const { return _bandData.begin(); }
		
		iterator end() { return _bandData.end(); }
		const_iterator end() const { return _bandData.end(); }
		
	private:
		std::vector<size_t> _dataDescToBand;
		std::vector<BandData> _bandData;
};

#endif
