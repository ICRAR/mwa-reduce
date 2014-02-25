#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include <fstream>
#include <string>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "../polarizationenum.h"
#include "../uvector.h"

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
	
	static Handle Partition(const string& msPath, size_t channelParts, class MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, const std::set<PolarizationEnum>& polsOut);
	
	class Handle {
	public:
		friend class PartitionedMS;
		
		Handle(const Handle& handle) : _metaFile(handle._metaFile), _msPath(handle._msPath), _channelParts(handle._channelParts), _polarizations(handle._polarizations), _referenceCount(handle._referenceCount)
		{
			++(*_referenceCount);
		}
		~Handle() { decrease(); }
		void operator=(const Handle& handle)
		{
			if(handle._referenceCount != _referenceCount)
			{
				decrease();
				_metaFile = handle._metaFile;
				_msPath = handle._msPath;
				_channelParts = handle._channelParts;
				_polarizations = handle._polarizations;
				_referenceCount = handle._referenceCount;
				++(*_referenceCount);
			}
		}
	private:
		void decrease();
		Handle(const std::string& metaFile, const std::string& msPath, size_t channelParts, const std::set<PolarizationEnum>& polarizations) : _metaFile(metaFile), _msPath(msPath), _channelParts(channelParts), _polarizations(polarizations), _referenceCount(new size_t(1)) { }
		std::string _metaFile, _msPath;
		size_t _channelParts;
		std::set<PolarizationEnum> _polarizations;
		size_t *_referenceCount;
	}; 
private:
	casa::MeasurementSet _ms;
	std::ifstream _metaFile, _weightFile, _dataFile;
	char *_modelFileMap;
	size_t _currentRow;
	bool _readPtrIsOk, _metaPtrIsOk, _weightPtrIsOk;
	ao::uvector<float> _weightBuffer;
	ao::uvector<std::complex<float>> _modelBuffer;
	
	struct MetaHeader
	{
		uint64_t selectedRowCount;
		uint32_t filenameLength;
		uint32_t fill;
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
	
	static std::string getPartPrefix(const std::string& msPath, size_t partIndex, PolarizationEnum pol);
	static std::string getMetaFilename(const std::string& msPath);
};

#endif
