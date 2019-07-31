#include "multibanddata.h"

MultiBandData::MultiBandData(casacore::MSSpectralWindow& spwTable, casacore::MSDataDescription& dataDescTable) :
_dataDescToBand(dataDescTable.nrow()),
_bandData(spwTable.nrow())
{
	for(size_t spw=0; spw!=_bandData.size(); ++spw)
	{
		_bandData[spw] = BandData(spwTable, spw);
	}
	
	casacore::ScalarColumn<int> spwColumn(dataDescTable, casacore::MSDataDescription::columnName(casacore::MSDataDescriptionEnums::SPECTRAL_WINDOW_ID));
	for(size_t id=0; id!=_dataDescToBand.size(); ++id)
		_dataDescToBand[id] = spwColumn(id);
}

MultiBandData::MultiBandData(const MultiBandData& source, size_t startChannel, size_t endChannel) :
	_dataDescToBand(source._dataDescToBand),
	_bandData(source.BandCount())
{
	for(size_t spw=0; spw!=source.BandCount(); ++spw)
		_bandData[spw] = BandData(source._bandData[spw], startChannel, endChannel);
}

std::set<size_t> MultiBandData::GetUsedDataDescIds(casacore::MeasurementSet& mainTable) const
{
	// If there is only one band, we assume it is used so as to avoid
	// scanning through the measurement set
	std::set<size_t> usedDataDescIds;
	if(_bandData.size() == 1)
		usedDataDescIds.insert(0);
	else
	{
		casacore::ScalarColumn<int> dataDescIdCol(mainTable, casacore::MeasurementSet::columnName(casacore::MSMainEnums::DATA_DESC_ID));
		for(size_t row=0; row!=mainTable.nrow(); ++row)
		{
			size_t dataDescId = dataDescIdCol(row);
			if(usedDataDescIds.find(dataDescId) == usedDataDescIds.end())
				usedDataDescIds.insert(dataDescId);
		}
	}
	return usedDataDescIds;
}
