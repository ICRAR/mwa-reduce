#ifndef MULTIBANDDATA_H
#define MULTIBANDDATA_H

#include <stdexcept>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include "banddata.h"

class MultiBandData
{
	public:
		MultiBandData()
		{
		}
		
		MultiBandData(casacore::MSSpectralWindow& spwTable, casacore::MSDataDescription& dataDescTable) :
		_dataDescToBand(dataDescTable.nrow()),
		_bandData(spwTable.nrow())
		{
			for(size_t spw=0; spw!=_bandData.size(); ++spw)
			{
				_bandData[spw] = BandData(spwTable, spw);
			}
			
			casacore::ROScalarColumn<int> spwColumn(	dataDescTable, casacore::MSDataDescription::columnName(casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));
			for(size_t id=0; id!=_dataDescToBand.size(); ++id)
				_dataDescToBand[id] = spwColumn(id);
		}
		
		MultiBandData(const MultiBandData& source) :
			_dataDescToBand(source._dataDescToBand),
			_bandData(source._bandData)
		{
		}
		
		MultiBandData(const MultiBandData& source, size_t startChannel, size_t endChannel) :
			_dataDescToBand(source._dataDescToBand)
		{
			_bandData.resize(source.BandCount());
			for(size_t spw=0; spw!=source.BandCount(); ++spw)
				_bandData[spw] = BandData(source._bandData[spw], startChannel, endChannel);
		}
		
		void operator=(const MultiBandData& source)
		{
			_dataDescToBand = source._dataDescToBand;
			_bandData = source._bandData;
		}
		
		const BandData& operator[](size_t dataDescID) const
		{
			return _bandData[_dataDescToBand[dataDescID]];
		}
		
		const BandData& FirstBand() const
		{
			return _bandData.front();
		}
		
		size_t BandCount() const
		{
			return _bandData.size();
		}
		
		size_t DataDescCount() const
		{
			return _dataDescToBand.size();
		}
		
		double LowestFrequency() const
		{
			double freq = _bandData[0].LowestFrequency();
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::min(freq, _bandData[i].LowestFrequency());
			return freq;
		}
		
		double CentreFrequency() const
		{
			return (LowestFrequency() + HighestFrequency()) * 0.5;
		}
		
		double HighestFrequency() const
		{
			double freq = _bandData[0].HighestFrequency();
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::max(freq, _bandData[i].HighestFrequency());
			return freq;
		}
		
		double Bandwidth() const
		{
			return HighestFrequency() - LowestFrequency();
		}
		
		double BandStart() const
		{
			double freq = std::min(_bandData[0].BandStart(), _bandData[0].BandEnd());
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::min(freq, std::min(_bandData[i].BandStart(), _bandData[i].BandEnd()));
			return freq;
		}
		
		double BandEnd() const
		{
			double freq = std::max(_bandData[0].BandStart(), _bandData[0].BandEnd());
			for(size_t i=0; i!=_bandData.size(); ++i)
				freq = std::max(freq, std::max(_bandData[i].BandStart(), _bandData[i].BandEnd()));
			return freq;
		}
		
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
	private:
		std::vector<size_t> _dataDescToBand;
		std::vector<BandData> _bandData;
};

#endif
