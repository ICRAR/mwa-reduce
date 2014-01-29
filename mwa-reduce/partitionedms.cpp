#include "partitionedms.h"
#include "multibanddata.h"

#include <cstdio>
#include <fstream>
#include <sstream>

#include <tables/Tables/ArrColDesc.h>

#define REDUNDANT_VALIDATION 1

ContiguousMS::ContiguousMS(const string& msPath, MSSelection selection, PolarizationEnum polOut, bool includeModel) :
	_timestep(0),
	_time(0.0),
	_dataDescId(0),
	_isModelColumnPrepared(false),
	_selection(selection),
	_polOut(polOut),
	_ms(msPath),
	_antenna1Column(_ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1)),
	_antenna2Column(_ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2)),
	_fieldIdColumn(_ms, casa::MS::columnName(casa::MSMainEnums::FIELD_ID)),
	_dataDescIdColumn(_ms, casa::MS::columnName(casa::MSMainEnums::DATA_DESC_ID)),
	_timeColumn(_ms, casa::MS::columnName(casa::MSMainEnums::TIME)),
	_uvwColumn(_ms, casa::MS::columnName(casa::MSMainEnums::UVW)),
	_weightColumn(_ms, casa::MS::columnName(casa::MSMainEnums::WEIGHT_SPECTRUM)),
	_dataColumn(_ms, casa::MS::columnName(casa::MSMainEnums::DATA)),
	_flagColumn(_ms, casa::MS::columnName(casa::MSMainEnums::FLAG))
{
	std::cout << "Opening " << msPath << " with contiguous MS reader.\n";
	const casa::IPosition shape(_dataColumn.shape(0));
	_polarizationCount = shape[0];
	_dataArray = casa::Array<std::complex<float>>(shape);
	_weightArray = casa::Array<float>(shape);
	_flagArray = casa::Array<bool>(shape);
	_bandData = MultiBandData(_ms.spectralWindow(), _ms.dataDescription());
	
	bool isWeightDefined = _weightColumn.isDefined(0);
	_msHasWeights = false;
	if(isWeightDefined)
	{
		casa::IPosition modelShape = _weightColumn.shape(0);
		_msHasWeights = (modelShape == shape);
	}
	if(!_msHasWeights)
	{
		_weightArray.set(1);
		std::cout << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM column; all visibilities are assumed to have equal weight.\n";
	}
	
	getRowRange(_ms, selection, _startRow, _endRow);
	Reset();
}

void ContiguousMS::Reset()
{
	_row = _startRow - 1;
	NextRow();
}

bool ContiguousMS::NextRow()
{
	int fieldId, a1, a2;
	do {
		++_row;
		if(_row >= _endRow)
			return false;
		
		fieldId = _fieldIdColumn(_row);
		a1 = _antenna1Column(_row);
		a2 = _antenna2Column(_row);
		if(_time != _timeColumn(_row))
		{
			++_timestep;
			_time = _timeColumn(_row);
		}
	} while(!_selection.IsSelected(fieldId, _timestep, a1, a2));
	
	_isMetaRead = false;
	_isDataRead = false;
	_isWeightRead = false;
	_isModelRead = false;
	
	return true;
}

void ContiguousMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId)
{
	readMeta();
	
	casa::Vector<double> uvwArray = _uvwColumn(_row);
	u = uvwArray(0);
	v = uvwArray(1);
	w = uvwArray(2);
	dataDescId = _dataDescId;
}

void ContiguousMS::ReadData(std::complex<float>* buffer)
{
	readMeta();
	readData();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeightedData(buffer,  startChannel, endChannel, _polarizationCount, _dataArray, _weightArray, _flagArray, _polOut);
}

void ContiguousMS::prepareModelColumn()
{
	if(_ms.isColumn(casa::MSMainEnums::MODEL_DATA))
	{
		_modelColumn.reset(new casa::ArrayColumn<std::complex<float> >(_ms, _ms.columnName(casa::MSMainEnums::MODEL_DATA)));
		casa::IPosition dataShape = _dataColumn.shape(0);
		bool isDefined = _modelColumn->isDefined(0);
		bool isSameShape = false;
		if(isDefined)
		{
			casa::IPosition modelShape = _modelColumn->shape(0);
			isSameShape = modelShape == dataShape;
		}
		if(!isDefined || !isSameShape)
		{
			std::cout << "WARNING: Your model column does not have the same shape as your data column: resetting MODEL column.\n";
			casa::Array<casa::Complex> zeroArray(dataShape);
			for(casa::Array<casa::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
				*i = std::complex<float>(0.0, 0.0);
			for(size_t row=0; row!=_ms.nrow(); ++row)
				_modelColumn->put(row, zeroArray);
		}
	}
	else { //if(!_ms.isColumn(casa::MSMainEnums::MODEL_DATA))
		std::cout << "Adding model data column... " << std::flush;
		casa::IPosition shape = _dataColumn.shape(0);
		casa::ArrayColumnDesc<casa::Complex> modelColumnDesc(_ms.columnName(casa::MSMainEnums::MODEL_DATA), shape);
		try {
			_ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
		} catch(std::exception& e)
		{
			_ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
		}
		
		casa::Array<casa::Complex> zeroArray(shape);
		for(casa::Array<casa::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
			*i = std::complex<float>(0.0, 0.0);
		
		_modelColumn.reset(new casa::ArrayColumn<std::complex<float> >(_ms, _ms.columnName(casa::MSMainEnums::MODEL_DATA)));
		
		for(size_t row=0; row!=_ms.nrow(); ++row)
			_modelColumn->put(row, zeroArray);
		
		std::cout << "DONE\n";
	}		
	const casa::IPosition shape(_modelColumn->shape(0));
	_modelArray = casa::Array<std::complex<float>>(shape);
	_isModelColumnPrepared = true;
}

void ContiguousMS::ReadModel(std::complex<float>* buffer)
{
	if(!_isModelColumnPrepared)
		prepareModelColumn();
	
	readMeta();
	readModel();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeightedData(buffer,  startChannel, endChannel, _polarizationCount, _modelArray, _weightArray, _flagArray, _polOut);
}

void ContiguousMS::WriteModel(size_t rowId, std::complex<float>* buffer)
{
	if(!_isModelColumnPrepared)
		prepareModelColumn();
	
	size_t dataDescId = _dataDescIdColumn(rowId);
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[dataDescId].ChannelCount();
	}
	
	_modelColumn->get(rowId, _modelArray);
	reverseCopyData(_modelArray, startChannel, endChannel, _polarizationCount, buffer, _polOut);
	_modelColumn->put(rowId, _modelArray);
}

void ContiguousMS::ReadWeights(std::complex<float>* buffer)
{
	readMeta();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeights(buffer,  startChannel, endChannel, _polarizationCount, _dataArray, _weightArray, _flagArray, _polOut);
}

void ContiguousMS::ReadWeights(float* buffer)
{
	readMeta();
	readWeights();
	size_t startChannel, endChannel;
	if(_selection.HasChannelRange())
	{
		startChannel = _selection.ChannelRangeStart();
		endChannel = _selection.ChannelRangeEnd();
	}
	else {
		startChannel = 0;
		endChannel = _bandData[_dataDescId].ChannelCount();
	}
	copyWeights(buffer,  startChannel, endChannel, _polarizationCount, _dataArray, _weightArray, _flagArray, _polOut);
}

PartitionedMS::PartitionedMS(const Handle& handle, size_t partIndex) :
	_metaFile(handle._metaFile),
	_currentRow(0),
	_readPtrIsOnData(true),
	_readPtrIsOnModel(false),
	_metaPtrIsOk(true),
	_weightPtrIsOk(true)
{
	_metaFile.read(reinterpret_cast<char*>(&_metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(_metaHeader.filenameLength+1, char(0));
	_metaFile.read(msPath.data(), _metaHeader.filenameLength);
	std::cout << "Opening reordered part " << partIndex << " for " << msPath.data() << '\n';
	_ms = casa::MeasurementSet(msPath.data());
	
	_dataFile.open(getPartPrefix(msPath.data(), partIndex)+".tmp", std::ios::in | std::ios::out);
	if(_dataFile.bad())
		throw std::runtime_error("Error opening temporary data file");
	_dataFile.read(reinterpret_cast<char*>(&_partHeader), sizeof(PartHeader));
	std::cout << "Channels " << _partHeader.channelStart << "-" << _partHeader.channelStart+_partHeader.channelCount << " are in part.\n";
	
	_weightFile.open(getPartPrefix(msPath.data(), partIndex)+"-w.tmp", std::ios::in | std::ios::out);
	if(_weightFile.bad())
		throw std::runtime_error("Error opening temporary data file");
	_weightBuffer.resize(_partHeader.channelCount);
}

void PartitionedMS::Reset()
{
	_currentRow = 0;
	_metaFile.seekg(sizeof(MetaHeader) + _metaHeader.filenameLength, std::ios::beg);
	_dataFile.seekg(sizeof(PartHeader), std::ios::beg);
	_weightFile.seekg(0, std::ios::beg);
	_readPtrIsOnData = true;
	_readPtrIsOnModel = false;
	_metaPtrIsOk = true;
	_weightPtrIsOk = true;
}

bool PartitionedMS::NextRow()
{
	++_currentRow;
	if(_currentRow >= _metaHeader.selectedRowCount)
		return false;
	if(_readPtrIsOnData)
	{
		if(_partHeader.hasModel)
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>) * 2, std::ios::cur);
		else
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
	else if(_readPtrIsOnModel)
	{
		if(_partHeader.hasModel)
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
	_readPtrIsOnData = true;
	_readPtrIsOnModel = false;
	
	if(_metaPtrIsOk)
		_metaFile.seekg(sizeof(MetaRecord), std::ios::cur);
	else
		_metaPtrIsOk = true;
	
	if(_weightPtrIsOk && _partHeader.hasWeights)
		_weightFile.seekg(_partHeader.channelCount * sizeof(float), std::ios::cur);
	_weightPtrIsOk = true;
	return true;
}

void PartitionedMS::ReadMeta(double& u, double& v, double& w, size_t& dataDescId)
{
	if(!_metaPtrIsOk)
		_metaFile.seekg(-sizeof(MetaRecord), std::ios::cur);
	_metaPtrIsOk = false;
	
	MetaRecord record;
	_metaFile.read(reinterpret_cast<char*>(&record), sizeof(MetaRecord));
	u = record.u;
	v = record.v;
	w = record.w;
	dataDescId = record.dataDescId;
}

void PartitionedMS::ReadData(std::complex<float>* buffer)
{
	if(!_readPtrIsOnData)
	{
		if(_readPtrIsOnModel)
			_dataFile.seekg(-_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
		else
			_dataFile.seekg(-2 * _partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
#ifdef REDUNDANT_VALIDATION
	size_t pos = size_t(_dataFile.tellg()) - sizeof(PartHeader);
	size_t fact = _partHeader.hasModel ? 2 : 1;
	if(pos != fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>))
	{
		std::ostringstream s;
		s << "Not on right pos: " << pos << " instead of " << fact * _currentRow * _partHeader.channelCount * sizeof(std::complex<float>) <<
			" (row " << (pos / (fact * _partHeader.channelCount * sizeof(std::complex<float>))) << " instead of " << _currentRow << ")";
		throw std::runtime_error(s.str());
	}
#endif
	_dataFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(std::complex<float>));
	_readPtrIsOnData = false;
	_readPtrIsOnModel = true;
}

void PartitionedMS::ReadModel(std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	if(!_readPtrIsOnModel)
	{
		if(_readPtrIsOnData)
			_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
		else
			_dataFile.seekg(-_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	}
	_dataFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(std::complex<float>));
	_readPtrIsOnData = false;
	_readPtrIsOnModel = false;
}

void PartitionedMS::WriteModel(size_t rowId, std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	
	_weightFile.seekg(_partHeader.channelCount * sizeof(float) * rowId, std::ios::beg);
	_weightFile.read(reinterpret_cast<char*>(_weightBuffer.data()), _partHeader.channelCount * sizeof(float));
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
	{
		buffer[i] *= _weightBuffer[i];
	}
	
	size_t pos = sizeof(PartHeader) + _partHeader.channelCount * sizeof(std::complex<float>) * (2 * rowId + 1);
	_dataFile.seekp(pos, std::ios::beg);
	_dataFile.write(reinterpret_cast<const char*>(buffer), _partHeader.channelCount * sizeof(std::complex<float>));
}

void PartitionedMS::ReadWeights(std::complex<float>* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	float* displacedBuffer = reinterpret_cast<float*>(buffer)+_partHeader.channelCount;
	_weightFile.read(reinterpret_cast<char*>(displacedBuffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
	copyRealToComplex(buffer, displacedBuffer, _partHeader.channelCount);
}

void PartitionedMS::ReadWeights(float* buffer)
{
	if(!_weightPtrIsOk)
		_weightFile.seekg(-_partHeader.channelCount * sizeof(float), std::ios::cur);
	_weightFile.read(reinterpret_cast<char*>(buffer), _partHeader.channelCount * sizeof(float));
	_weightPtrIsOk = false;
}

std::string PartitionedMS::getPartPrefix(const std::string& msPath, size_t partIndex)
{
	std::string prefix(msPath);
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	
	std::ostringstream partPrefix;
	partPrefix << prefix << "-part";
	if(partIndex < 1000) partPrefix << '0';
	if(partIndex < 100) partPrefix << '0';
	if(partIndex < 10) partPrefix << '0';
	partPrefix << partIndex;
	return partPrefix.str();
}

string PartitionedMS::getMetaFilename(const string& msPath)
{
	std::string prefix(msPath);
	while(!prefix.empty() && *prefix.rbegin() == '/')
		prefix.resize(prefix.size()-1);
	return prefix + "-parted-meta.tmp";
}

/*
 * When partitioned:
 * One global file stores:
 * - Metadata:
 *   * Number of selected rows
 *   * Filename length + string
 *   * [ UVW, dataDescId ]
 * The binary parts store the following information:
 * - Number of channels
 * - Start channel in MS
 * - Total weight in part
 * - Data    (single polarization, as requested)
 * - Weights (single, only needed when imaging PSF)
 * - Model, optionally
 */
PartitionedMS::Handle PartitionedMS::Partition(const string& msPath, size_t channelParts, MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, PolarizationEnum polOut)
{
	casa::MeasurementSet ms(msPath);

	std::vector<std::ofstream*> dataFiles(channelParts), weightFiles(channelParts);
	for(size_t part=0; part!=channelParts; ++part)
	{
		std::string partPrefix = getPartPrefix(msPath, part);
		dataFiles[part] = new std::ofstream(partPrefix + ".tmp");
		if(includeWeights)
			weightFiles[part] = new std::ofstream(partPrefix + "-w.tmp");
		dataFiles[part]->seekp(sizeof(PartHeader), std::ios::beg);
	}
	MultiBandData band(ms.spectralWindow(), ms.dataDescription());
	casa::ROScalarColumn<int> antenna1Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> antenna2Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROScalarColumn<int> fieldIdColumn(ms, casa::MS::columnName(casa::MSMainEnums::FIELD_ID));
	casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
	casa::ROArrayColumn<double> uvwColumn(ms, casa::MS::columnName(casa::MSMainEnums::UVW));
	casa::ROArrayColumn<float> weightColumn(ms, casa::MS::columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, dataColumnName);
	casa::ROArrayColumn<bool> flagColumn(ms, casa::MS::columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<int> dataDescIdColumn(ms, ms.columnName(casa::MSMainEnums::DATA_DESC_ID));
	
	const casa::IPosition shape(dataColumn.shape(0));
	const size_t polarizationCount = shape[0];
	size_t channelCount, channelStart;
	if(selection.HasChannelRange())
	{
		channelCount = selection.ChannelRangeEnd() - selection.ChannelRangeStart();
		channelStart = selection.ChannelRangeStart();
	}
	else {
		channelCount = shape[1];
		channelStart = 0;
	}
	
	size_t startRow, endRow;
	getRowRange(ms, selection, startRow, endRow);
	
	// Count selected rows
	uint64_t selectedRowCount = 0;
	size_t timestep = 0;
	double time = timeColumn(0);
	for(size_t row=startRow; row!=endRow; ++row)
	{
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row), fieldId = fieldIdColumn(row);
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		if(selection.IsSelected(fieldId, timestep, a1, a2))
			++selectedRowCount;
	}
	std::cout << "Reordering " << selectedRowCount << " selected rows into " << channelParts << " parts.\n";

	// Write header of meta file
	std::string metaFilename = getMetaFilename(msPath);
	std::ofstream metaFile(metaFilename);
	MetaHeader metaHeader;
	metaHeader.selectedRowCount = selectedRowCount;
	metaHeader.filenameLength = msPath.size();
	metaHeader.fill = 0;
	metaFile.write(reinterpret_cast<char*>(&metaHeader), sizeof(metaHeader));
	metaFile.write(msPath.c_str(), msPath.size());
	
	// Write actual data
	timestep = 0;
	time = timeColumn(0);
	
	std::vector<std::complex<float>> dataBuffer(polarizationCount * (1 + channelCount / channelParts));
	std::vector<float> weightBuffer(polarizationCount * (1 + channelCount / channelParts));
	
	casa::Array<std::complex<float>> dataArray(shape);
	casa::Array<float> weightArray(shape);
	casa::Array<bool> flagArray(shape);
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
			
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		if(selection.IsSelected(fieldId, timestep, a1, a2))
		{
			size_t dataDescId = dataDescIdColumn(row);
			casa::Vector<double> uvwArray = uvwColumn(row);
			MetaRecord meta;
			meta.u = uvwArray(0);
			meta.v = uvwArray(1);
			meta.w = uvwArray(2);
			meta.dataDescId = dataDescId;
			metaFile.write(reinterpret_cast<char*>(&meta), sizeof(MetaRecord));
			if(metaFile.bad())
				throw std::runtime_error("Error writing to temporary file");
				
			dataColumn.get(row, dataArray);
			weightColumn.get(row, weightArray);
			flagColumn.get(row, flagArray);
			
			for(size_t part=0; part!=channelParts; ++part)
			{
				size_t
					partStartCh = channelStart + channelCount*part/channelParts,
					partEndCh = channelStart + channelCount*(part+1)/channelParts;
				
				copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, polarizationCount, dataArray, weightArray, flagArray, polOut);
				dataFiles[part]->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
				
				if(includeModel)
				{
					// The initial model in the MS is never used, so we just have to reserve room for later
					dataFiles[part]->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
				}
				if(dataFiles[part]->bad())
					throw std::runtime_error("Error writing to temporary data file");
				
				if(includeWeights)
				{
					copyWeights(weightBuffer.data(), partStartCh, partEndCh, polarizationCount, dataArray, weightArray, flagArray, polOut);
					weightFiles[part]->write(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
					if(weightFiles[part]->bad())
						throw std::runtime_error("Error writing to temporary weights file");
				}
			}
		}
	}
	
	// Write header to parts
	for(size_t part=0; part!=channelParts; ++part)
	{
		PartHeader header;
		header.channelStart = channelStart + channelCount*part/channelParts,
		header.channelCount = (channelStart + channelCount*(part+1)/channelParts) - header.channelStart;
		header.hasModel = includeModel;
		header.hasWeights = includeWeights;
		dataFiles[part]->seekp(0, std::ios::beg);
		dataFiles[part]->write(reinterpret_cast<char*>(&header), sizeof(PartHeader));
		if(dataFiles[part]->bad())
			throw std::runtime_error("Error writing to temporary data file");
		
		delete dataFiles[part];
		if(includeWeights)
			delete weightFiles[part];
	}
	
	return Handle(metaFilename, msPath, channelParts);
}

void MSProvider::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(polOut == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			} else {
				dest[ch] = 0;
			}
			weightPtr += polCount-1;
			inPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] += *inPtr * (*weightPtr);
			}
			++weightPtr;
			++inPtr;
			++flagPtr;
		}
	} /*else if(Polarization() == Polarization::XY || Polarization() == Polarization::YX)
	{
		// Step to XY:
		++weightPtr;
		++inPtr;
		++flagPtr;
		const bool flipSign = Polarization() == Polarization::YX;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			// Step to YX:
			++weightPtr;
			++inPtr;
			++flagPtr;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				//dest[ch] += *inPtr * (*weightPtr) * rowWeight;
				//_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 3;
			inPtr += 3;
			flagPtr += 3;
			//if(flipSign)
			//	dest[ch].imag(dest[ch].imag() * -1.0);
		}
	}*/ else {
		int polIndex = Polarization::TypeToIndex(polOut, polCount);
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			}
			else {
				dest[ch] = 0;
			}
			weightPtr += polCount;
			inPtr += polCount;
			flagPtr += polCount;
		}
	}
}

template<typename NumType>
void MSProvider::copyWeights(NumType* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(polOut == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polCount-1;
			weightPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] += *weightPtr;
			++inPtr;
			++weightPtr;
			++flagPtr;
		}
	} /*else if(Polarization() == Polarization::XY || Polarization() == Polarization::YX)
	{
		// Step to XY:
		inPtr++;
		weightPtr++;
		flagPtr++;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			// Step to YX:
			inPtr++;
			weightPtr++;
			flagPtr++;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				//dest[ch] += (*weightPtr) * rowWeight;
				//_totalWeight += (*weightPtr) * rowWeight;
			}
			inPtr += 3;
			weightPtr += 3;
			flagPtr += 3;
		}
	} */ else {
		int polIndex = Polarization::TypeToIndex(polOut, polCount);
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polCount;
			weightPtr += polCount;
			flagPtr += polCount;
		}
	}
}

template
void MSProvider::copyWeights<float>(float* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

template
void MSProvider::copyWeights<std::complex<float>>(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

void MSProvider::reverseCopyData(casa::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, size_t polCount, const std::complex<float>* source, PolarizationEnum polSource)
{
	const size_t selectedChannelCount = endChannel - startChannel;
	casa::Array<std::complex<float>>::contiter dataIter = dest.cbegin() + startChannel * polCount;
	
	if(polSource == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*dataIter = source[ch];
				*(dataIter + (polCount-1)) = source[ch];
			}
			dataIter += polCount;
		}
	} else {
		int polIndex = Polarization::TypeToIndex(polSource, polCount);
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*(dataIter+polIndex) = source[ch];
			}
			dataIter += polCount;
		}
	}
}

void MSProvider::getRowRange(casa::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow)
{
	startRow = 0;
	endRow = ms.nrow();
	if(selection.HasInterval())
	{
		std::cout << "Determining first and last row index... " << std::flush;
		casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
		double time = timeColumn(0);
		size_t timestepIndex = 0;
		for(size_t row = 0; row!=ms.nrow(); ++row)
		{
			if(time != timeColumn(row))
			{
				++timestepIndex;
				if(timestepIndex == selection.IntervalStart())
					startRow = row;
				if(timestepIndex == selection.IntervalEnd())
				{
					endRow = row;
					break;
				}
				time = timeColumn(row);
			}
		}
		std::cout << "DONE (" << startRow << '-' << endRow << ")\n";
	}
}

void PartitionedMS::Handle::decrease()
{
	--(*_referenceCount);
	if((*_referenceCount) == 0)
	{
		std::cout << "Cleaning up temporary files...\n";
		for(size_t part=0; part!=_channelParts; ++part)
		{
			std::string prefix = getPartPrefix(_msPath, part);
			std::remove((prefix + ".tmp").c_str());
			std::remove((prefix + "-w.tmp").c_str());
		}
		std::remove(_metaFile.c_str());
		delete _referenceCount;
	}
}
