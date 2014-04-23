#include "partitionedms.h"

#include "../multibanddata.h"
#include "../progressbar.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <memory>
// #define REDUNDANT_VALIDATION 1

PartitionedMS::PartitionedMS(const Handle& handle, size_t partIndex, PolarizationEnum polarization) :
	_metaFile(handle._data->_metaFile),
	_modelFileMap(0),
	_currentRow(0),
	_readPtrIsOk(true),
	_metaPtrIsOk(true),
	_weightPtrIsOk(true)
{
	_metaFile.read(reinterpret_cast<char*>(&_metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(_metaHeader.filenameLength+1, char(0));
	_metaFile.read(msPath.data(), _metaHeader.filenameLength);
	std::cout << "Opening reordered part " << partIndex << " for " << msPath.data() << '\n';
	_ms = casa::MeasurementSet(msPath.data());
	
	_dataFile.open(getPartPrefix(msPath.data(), partIndex, polarization)+".tmp", std::ios::in);
	if(_dataFile.bad())
		throw std::runtime_error("Error opening temporary data file");
	_dataFile.read(reinterpret_cast<char*>(&_partHeader), sizeof(PartHeader));
	
	if(_partHeader.hasModel)
	{
		_fd = open((getPartPrefix(msPath.data(), partIndex, polarization)+"-m.tmp").c_str(), O_RDWR);
		if(_fd == -1)
			throw std::runtime_error("Error opening temporary data file");
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		_modelFileMap = reinterpret_cast<char*>( mmap(NULL, length, PROT_WRITE | PROT_WRITE, MAP_SHARED, _fd, 0) );
		if(_modelFileMap == MAP_FAILED)
		{
			_modelFileMap = 0;
			throw std::runtime_error("Error creating memory map to temporary model file: mmap() return MAP_FAILED.");
		}
	}
	
	_weightFile.open(getPartPrefix(msPath.data(), partIndex, polarization)+"-w.tmp", std::ios::in);
	if(_weightFile.bad())
		throw std::runtime_error("Error opening temporary data file");
	_weightBuffer.resize(_partHeader.channelCount);
	_modelBuffer.resize(_partHeader.channelCount);
}

PartitionedMS::~PartitionedMS()
{
	if(_modelFileMap!=0)
	{
		size_t length = _partHeader.channelCount * _metaHeader.selectedRowCount * sizeof(std::complex<float>);
		munmap(_modelFileMap, length);
	}
	if(_partHeader.hasModel)
		close(_fd);
}

void PartitionedMS::Reset()
{
	_currentRow = 0;
	_metaFile.seekg(sizeof(MetaHeader) + _metaHeader.filenameLength, std::ios::beg);
	_dataFile.seekg(sizeof(PartHeader), std::ios::beg);
	_weightFile.seekg(0, std::ios::beg);
	_readPtrIsOk = true;
	_metaPtrIsOk = true;
	_weightPtrIsOk = true;
}

bool PartitionedMS::NextRow()
{
	++_currentRow;
	if(_currentRow >= _metaHeader.selectedRowCount)
		return false;
	
	if(_readPtrIsOk)
		_dataFile.seekg(_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
	else
		_readPtrIsOk = true;
	
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
	if(!_readPtrIsOk)
	{
		_dataFile.seekg(-_partHeader.channelCount * sizeof(std::complex<float>), std::ios::cur);
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
	_readPtrIsOk = false;
}

void PartitionedMS::ReadModel(std::complex<float>* buffer)
{
#ifdef REDUNDANT_VALIDATION
	if(!_partHeader.hasModel)
		throw std::runtime_error("Partitioned MS initialized without model");
#endif
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	memcpy(reinterpret_cast<char*>(buffer), _modelFileMap + rowLength*_currentRow, rowLength);
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
		buffer[i] *= _weightBuffer[i];
	
	size_t rowLength = _partHeader.channelCount * sizeof(std::complex<float>);
	std::complex<float>* modelWritePtr = reinterpret_cast<std::complex<float>*>(_modelFileMap + rowLength*rowId);
	
	// In case the value was not sampled in this pass, it will be set to infinite and should not overwrite the current
	// value in the set.
	for(size_t i=0; i!=_partHeader.channelCount; ++i)
	{
		if(std::isfinite(buffer[i].real()))
			modelWritePtr[i] = buffer[i];
	}
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

std::string PartitionedMS::getPartPrefix(const std::string& msPath, size_t partIndex, PolarizationEnum pol)
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
	partPrefix << "-";
	partPrefix << Polarization::TypeToShortString(pol);
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
PartitionedMS::Handle PartitionedMS::Partition(const string& msPath, size_t channelParts, MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, const std::set<PolarizationEnum>& polsOut)
{
	casa::MeasurementSet ms(msPath);

	std::vector<std::ofstream*> dataFiles(channelParts*polsOut.size()), weightFiles(channelParts*polsOut.size());
	size_t fileIndex = 0;
	for(size_t part=0; part!=channelParts; ++part)
	{
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			std::string partPrefix = getPartPrefix(msPath, part, *p);
			dataFiles[fileIndex] = new std::ofstream(partPrefix + ".tmp");
			if(includeWeights)
				weightFiles[fileIndex] = new std::ofstream(partPrefix + "-w.tmp");
			dataFiles[fileIndex]->seekp(sizeof(PartHeader), std::ios::beg);
			++fileIndex;
		}
	}
	std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(ms);
	
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
	
	bool isWeightDefined = weightColumn.isDefined(0);
	bool msHasWeights = false;
	casa::Array<float> weightArray(shape);
	if(isWeightDefined)
	{
		casa::IPosition modelShape = weightColumn.shape(0);
		msHasWeights = (modelShape == shape);
	}
	if(!msHasWeights)
	{
		weightArray.set(1);
		std::cout << "WARNING: This measurement set has no or an invalid WEIGHT_SPECTRUM column; all visibilities are assumed to have equal weight.\n";
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
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
		casa::Vector<double> uvw = uvwColumn(row);
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			++selectedRowCount;
	}
	std::cout << "Reordering " << selectedRowCount << " selected rows into " << channelParts << " x " << polsOut.size() << " parts.\n";

	// Write header of meta file
	std::string metaFilename = getMetaFilename(msPath);
	std::ofstream metaFile(metaFilename);
	MetaHeader metaHeader;
	memset(&metaHeader, 0, sizeof(MetaHeader));
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
	casa::Array<bool> flagArray(shape);
	ProgressBar progress1("Reordering");
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		progress1.SetProgress(row, ms.nrow());
		const int
			a1 = antenna1Column(row), a2 = antenna2Column(row),
			fieldId = fieldIdColumn(row);
			
		if(time != timeColumn(row))
		{
			++timestep;
			time = timeColumn(row);
		}
		casa::Vector<double> uvwArray = uvwColumn(row);
		if(selection.IsSelected(fieldId, timestep, a1, a2, uvwArray))
		{
			size_t dataDescId = dataDescIdColumn(row);
			MetaRecord meta;
			memset(&meta, 0, sizeof(MetaRecord));
			meta.u = uvwArray(0);
			meta.v = uvwArray(1);
			meta.w = uvwArray(2);
			meta.dataDescId = dataDescId;
			metaFile.write(reinterpret_cast<char*>(&meta), sizeof(MetaRecord));
			if(metaFile.bad())
				throw std::runtime_error("Error writing to temporary file");
				
			dataColumn.get(row, dataArray);
			if(msHasWeights)
				weightColumn.get(row, weightArray);
			flagColumn.get(row, flagArray);
			
			fileIndex = 0;
			for(size_t part=0; part!=channelParts; ++part)
			{
				size_t
					partStartCh = channelStart + channelCount*part/channelParts,
					partEndCh = channelStart + channelCount*(part+1)/channelParts;
				
				for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
				{
					copyWeightedData(dataBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightArray, flagArray, *p);
					dataFiles[fileIndex]->write(reinterpret_cast<char*>(dataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
					if(dataFiles[fileIndex]->bad())
						throw std::runtime_error("Error writing to temporary data file");
					
					if(includeWeights)
					{
						copyWeights(weightBuffer.data(), partStartCh, partEndCh, msPolarizations, dataArray, weightArray, flagArray, *p);
						weightFiles[fileIndex]->write(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
						if(weightFiles[fileIndex]->bad())
							throw std::runtime_error("Error writing to temporary weights file");
					}
					++fileIndex;
				}
			}
		}
	}
	progress1.SetProgress(ms.nrow(), ms.nrow());
	
	// Write header to parts and write model files
	PartHeader header;
	memset(&header, 0, sizeof(PartHeader));
	header.hasModel = includeModel;
	header.hasWeights = includeWeights;
	fileIndex = 0;
	memset(reinterpret_cast<char*>(dataBuffer.data()), 0, sizeof(std::complex<float>)*(1 + channelCount / channelParts));
	std::unique_ptr<ProgressBar> progress2;
	if(includeModel)
		progress2.reset(new ProgressBar("Initializing model visibilities"));
	for(size_t part=0; part!=channelParts; ++part)
	{
		header.channelStart = channelStart + channelCount*part/channelParts,
		header.channelCount = (channelStart + channelCount*(part+1)/channelParts) - header.channelStart;
		for(std::set<PolarizationEnum>::const_iterator p=polsOut.begin(); p!=polsOut.end(); ++p)
		{
			dataFiles[fileIndex]->seekp(0, std::ios::beg);
			dataFiles[fileIndex]->write(reinterpret_cast<char*>(&header), sizeof(PartHeader));
			if(dataFiles[fileIndex]->bad())
				throw std::runtime_error("Error writing to temporary data file");
			
			delete dataFiles[fileIndex];
			if(includeWeights)
				delete weightFiles[fileIndex];
			++fileIndex;
			
			// If model is requested, fill model file with zeros
			if(includeModel)
			{
				std::string partPrefix = getPartPrefix(msPath, part, *p);
				std::ofstream modelFile(partPrefix + "-m.tmp");
				for(size_t i=0; i!=selectedRowCount; ++i)
				{
					modelFile.write(reinterpret_cast<char*>(dataBuffer.data()), header.channelCount * sizeof(std::complex<float>));
					progress2->SetProgress(part*selectedRowCount + i, channelParts*selectedRowCount);
				}
			}
		}
	}
	progress2.reset();
	
	return Handle(metaFilename, msPath, dataColumnName, channelParts, polsOut, selection);
}

void PartitionedMS::unpartition(const PartitionedMS::Handle& handle)
{
	const std::set<PolarizationEnum> pols = handle._data->_polarizations;
	std::ifstream metaFile(handle._data->_metaFile);
	MetaHeader metaHeader;
	metaFile.read(reinterpret_cast<char*>(&metaHeader), sizeof(MetaHeader));
	std::vector<char> msPath(metaHeader.filenameLength+1, char(0));
	metaFile.read(msPath.data(), metaHeader.filenameLength);
	
	std::ifstream firstDataFile(getPartPrefix(msPath.data(), 0, *pols.begin())+".tmp", std::ios::in);
	if(firstDataFile.bad())
		throw std::runtime_error("Error opening temporary data file");
	PartHeader firstPartHeader;
	firstDataFile.read(reinterpret_cast<char*>(&firstPartHeader), sizeof(PartHeader));
	
	if(firstPartHeader.hasModel)
	{
		const size_t channelParts = handle._data->_channelParts;
		std::vector<std::ifstream*> modelFiles(channelParts*pols.size()), weightFiles(channelParts*pols.size());
		size_t fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				std::string partPrefix = getPartPrefix(msPath.data(), part, *p);
				modelFiles[fileIndex] = new std::ifstream(partPrefix + "-m.tmp");
				if(firstPartHeader.hasWeights)
					weightFiles[fileIndex] = new std::ifstream(partPrefix + "-w.tmp");
				++fileIndex;
			}
		}
		
		casa::MeasurementSet ms(msPath.data(), casa::Table::Update);
		const std::vector<PolarizationEnum> msPolarizations = GetMSPolarizations(ms);
		initializeModelColumn(ms);
		casa::ROScalarColumn<int> antenna1Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> antenna2Column(ms, casa::MS::columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROScalarColumn<int> fieldIdColumn(ms, casa::MS::columnName(casa::MSMainEnums::FIELD_ID));
		casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
		casa::ROArrayColumn<casa::Complex> dataColumn(ms, handle._data->_dataColumnName);
		casa::ArrayColumn<casa::Complex> modelColumn(ms, casa::MS::columnName(casa::MSMainEnums::MODEL_DATA));
		casa::ROArrayColumn<double> uvwColumn(ms, casa::MS::columnName(casa::MSMainEnums::UVW));
		
		const casa::IPosition shape(dataColumn.shape(0));
		size_t channelCount, channelStart;
		if(handle._data->_selection.HasChannelRange())
		{
			channelCount = handle._data->_selection.ChannelRangeEnd() - handle._data->_selection.ChannelRangeStart();
			channelStart = handle._data->_selection.ChannelRangeStart();
		}
		else {
			channelCount = shape[1];
			channelStart = 0;
		}
		
		std::vector<std::complex<float>> modelDataBuffer(1 + channelCount / channelParts);
		std::vector<float> weightBuffer(1 + channelCount / channelParts);
		casa::Array<std::complex<float>> modelDataArray(shape);
	
		ProgressBar progress(std::string("Writing changed model back to ") + msPath.data());
		size_t timestep = 0;
		double time = timeColumn(0);
		for(size_t row=0; row!=ms.nrow(); ++row)
		{
			progress.SetProgress(row, ms.nrow());
			const int
				a1 = antenna1Column(row), a2 = antenna2Column(row),
				fieldId = fieldIdColumn(row);
			casa::Vector<double> uvw = uvwColumn(row);
				
			if(time != timeColumn(row))
			{
				++timestep;
				time = timeColumn(row);
			}
			if(handle._data->_selection.IsSelected(fieldId, timestep, a1, a2, uvw))
			{
				modelColumn.get(row, modelDataArray);
				size_t fileIndex = 0;
				for(size_t part=0; part!=channelParts; ++part)
				{
					size_t
						partStartCh = channelStart + channelCount*part/channelParts,
						partEndCh = channelStart + channelCount*(part+1)/channelParts;
					
					for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
					{
						modelFiles[fileIndex]->read(reinterpret_cast<char*>(modelDataBuffer.data()), (partEndCh - partStartCh) * sizeof(std::complex<float>));
						if(firstPartHeader.hasWeights)
						{
							weightFiles[fileIndex]->read(reinterpret_cast<char*>(weightBuffer.data()), (partEndCh - partStartCh) * sizeof(float));
							for(size_t i=0; i!=partEndCh - partStartCh; ++i)
								modelDataBuffer[i] /= weightBuffer[i];
						}
						if(modelFiles[fileIndex]->bad())
							throw std::runtime_error("Error writing to temporary data file");
						reverseCopyData(modelDataArray, partStartCh, partEndCh, msPolarizations, modelDataBuffer.data(), *p);
						
						++fileIndex;
					}
				}
				modelColumn.put(row, modelDataArray);
			}
		}
		progress.SetProgress(ms.nrow(),ms.nrow());
		
		fileIndex = 0;
		for(size_t part=0; part!=channelParts; ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=pols.begin(); p!=pols.end(); ++p)
			{
				delete modelFiles[fileIndex];
				if(firstPartHeader.hasWeights)
					delete weightFiles[fileIndex];
				++fileIndex;
			}
		}
	}
}

void PartitionedMS::Handle::decrease()
{
	--(_data->_referenceCount);
	if(_data->_referenceCount == 0)
	{
		PartitionedMS::unpartition(*this);
		
		std::cout << "Cleaning up temporary files...\n";
		
		//MetaHeader metaHeader;
		//std::ifstream metaFile(_metaFile);
		//metaFile.read(reinterpret_cast<char*>(&metaHeader), sizeof(MetaHeader));
		
		for(size_t part=0; part!=_data->_channelParts; ++part)
		{
			for(std::set<PolarizationEnum>::const_iterator p=_data->_polarizations.begin(); p!=_data->_polarizations.end(); ++p)
			{
				std::string prefix = getPartPrefix(_data->_msPath, part, *p);
				std::remove((prefix + ".tmp").c_str());
				std::remove((prefix + "-w.tmp").c_str());
				std::remove((prefix + "-m.tmp").c_str());
			}
		}
		std::remove(_data->_metaFile.c_str());
		delete _data;
	}
}
