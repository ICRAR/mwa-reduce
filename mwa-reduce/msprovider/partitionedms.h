#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include <fstream>
#include <string>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "../polarizationenum.h"
#include "../uvector.h"
#include "../msselection.h"

#include "msprovider.h"

class PartitionedMS : public MSProvider
{
public:
	class Handle;
	
	PartitionedMS(const Handle& handle, size_t partIndex, PolarizationEnum polarization);
	virtual ~PartitionedMS();
	
	virtual casa::MeasurementSet &MS() { return _ms; }
	
	virtual size_t RowId() const { return _currentRow; }
	
	virtual bool NextRow();
	
	virtual void Reset();
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId);
	
	virtual void ReadData(std::complex<float>* buffer);
	
	virtual void ReadModel(std::complex<float>* buffer);
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer);
	
	virtual void ReadWeights(float* buffer);
	
	virtual void ReadWeights(std::complex<float>* buffer);
	
	virtual void ReopenRW() { }
	
	virtual double StartTime() { return _metaHeader.startTime; }
	
	static Handle Partition(const string& msPath, size_t channelParts, class MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, const std::set<PolarizationEnum>& polsOut, const std::string& temporaryDirectory);
	
	class Handle {
	public:
		friend class PartitionedMS;
		
		Handle(const Handle& handle) : _data(handle._data)
		{
			++(_data->_referenceCount);
		}
		~Handle() { decrease(); }
		void operator=(const Handle& handle)
		{
			if(handle._data != _data)
			{
				decrease();
				_data = handle._data;
				++(_data->_referenceCount);
			}
		}
	private:
		struct HandleData
		{
			HandleData(const std::string& metaFile, const std::string& msPath, const string& dataColumnName, const std::string& temporaryDirectory, size_t channelParts, const std::set<PolarizationEnum>& polarizations, const MSSelection& selection) :
			_metaFile(metaFile), _msPath(msPath), _dataColumnName(dataColumnName), _temporaryDirectory(temporaryDirectory), _channelParts(channelParts), _polarizations(polarizations), _selection(selection), _referenceCount(1) { }
			
			std::string _metaFile, _msPath, _dataColumnName, _temporaryDirectory;
			size_t _channelParts;
			std::set<PolarizationEnum> _polarizations;
			MSSelection _selection;
			size_t _referenceCount;
		} *_data;
		
		void decrease();
		Handle(const std::string& metaFile, const std::string& msPath, const string& dataColumnName, const std::string& temporaryDirectory, size_t channelParts, const std::set<PolarizationEnum>& polarizations, const MSSelection& selection) :
			_data(new HandleData(metaFile, msPath, dataColumnName, temporaryDirectory, channelParts, polarizations, selection))
		{
		}
	}; 
private:
	static void unpartition(const Handle& handle);
	
	casa::MeasurementSet _ms;
	std::ifstream _metaFile, _weightFile, _dataFile;
	char *_modelFileMap;
	size_t _currentRow;
	bool _readPtrIsOk, _metaPtrIsOk, _weightPtrIsOk;
	ao::uvector<float> _weightBuffer;
	ao::uvector<std::complex<float>> _modelBuffer;
	int _fd;
	
	struct MetaHeader
	{
		uint64_t selectedRowCount;
		uint32_t filenameLength;
		double startTime;
	} _metaHeader;
	struct MetaRecord
	{
		double u, v, w;
		uint32_t dataDescId;
	};
	struct PartHeader
	{
		uint64_t channelCount;
		uint64_t channelStart;
		bool hasModel, hasWeights;
	} _partHeader;
	
	static std::string getPartPrefix(const std::string& msPath, size_t partIndex, PolarizationEnum pol, const std::string& tempDir);
	static std::string getMetaFilename(const std::string& msPath, const std::string& tempDir);
};

#endif
