#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "msselection.h"
#include "polarizationenum.h"
#include "multibanddata.h"

class MSProvider
{
public:
	virtual casa::MeasurementSet &MS() = 0;
	
	virtual size_t RowId() const = 0;
	
	virtual bool NextRow() = 0;
	
	virtual void Reset() = 0;
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId) = 0;
	
	virtual void ReadData(std::complex<float>* buffer) = 0;
	
	virtual void ReadModel(std::complex<float>* buffer) = 0;
	
	virtual void WriteModel(size_t rowId, const std::complex<float>* buffer) = 0;
	
	virtual void ReadWeights(float* buffer) = 0;
	
	virtual void ReadWeights(std::complex<float>* buffer) = 0;
	
	virtual double TotalWeight() = 0;
protected:
	static void copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut, double& totalWeight);
	
	template<typename NumType>
	static void copyWeights(NumType* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);
	
	static void reverseCopyData(casa::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, size_t polCount, const std::complex<float>* source, PolarizationEnum polSource, bool addData);
	
	static void getRowRange(casa::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow);
	
	static void copyRealToComplex(std::complex<float>* dest, const float* source, size_t n)
	{
		const float* end = source + n;
		while(source != end)
		{
			*dest = *source;
			++dest;
			++source;
		}
	}
	
	MSProvider() { }
private:
	MSProvider(const MSProvider&) { }
	void operator=(const MSProvider&) { }
};

class ContiguousMS : public MSProvider
{
public:
	ContiguousMS(const string& msPath, MSSelection selection, PolarizationEnum polOut, bool includeModel);
	
	virtual casa::MeasurementSet &MS() { return _ms; }
	
	virtual size_t RowId() const { return _row; }
	
	virtual bool NextRow();
	
	virtual void Reset();
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId);
	
	virtual void ReadData(std::complex<float>* buffer);
	
	virtual void ReadModel(std::complex<float>* buffer);
	
	virtual void WriteModel(size_t rowId, const std::complex<float>* buffer);
	
	virtual void ReadWeights(float* buffer);
	
	virtual void ReadWeights(std::complex<float>* buffer);
	
	virtual double TotalWeight();
private:
	size_t _row;
	size_t _timestep;
	double _time, _totalWeight;
	int _dataDescId;
	size_t _startRow, _endRow, _polarizationCount;
	MSSelection _selection;
	PolarizationEnum _polOut;
	casa::MeasurementSet _ms;
	MultiBandData _bandData;
	bool _msHasWeights;

	casa::ROScalarColumn<int> _antenna1Column, _antenna2Column, _fieldIdColumn, _dataDescIdColumn;
	casa::ROScalarColumn<double> _timeColumn;
	casa::ROArrayColumn<double> _uvwColumn;
	casa::ROArrayColumn<float> _weightColumn;
	casa::ROArrayColumn<casa::Complex> _dataColumn;
	casa::ROArrayColumn<bool> _flagColumn;
	std::unique_ptr<casa::ArrayColumn<casa::Complex>> _modelColumn;
	
	casa::Array<std::complex<float>> _dataArray, _modelArray;
	casa::Array<float> _weightArray;
	casa::Array<bool> _flagArray;
};

class PartitionedMS : public MSProvider
{
public:
	PartitionedMS(const string& metaFilename, size_t partIndex);
	
	virtual casa::MeasurementSet &MS() { return _ms; }
	
	virtual size_t RowId() const { return _currentRow; }
	
	virtual bool NextRow();
	
	virtual void Reset();
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId);
	
	virtual void ReadData(std::complex<float>* buffer);
	
	virtual void ReadModel(std::complex<float>* buffer);
	
	virtual void WriteModel(size_t rowId, const std::complex<float>* buffer);
	
	virtual void ReadWeights(float* buffer);
	
	virtual void ReadWeights(std::complex<float>* buffer);
	
	virtual double TotalWeight() { return _partHeader.totalWeight; }
	
	static std::string Partition(const string& msPath, size_t channelParts, MSSelection& selection, const string& dataColumnName, bool includeWeights, bool includeModel, PolarizationEnum polOut);
private:
	casa::MeasurementSet _ms;
	std::ifstream _metaFile, _weightFile;
	std::fstream _dataFile;
	size_t _currentRow;
	bool _readPtrIsOnData, _readPtrIsOnModel, _metaPtrIsOk, _weightPtrIsOk;
	
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
		double totalWeight;
		bool hasModel, hasWeights;
	} _partHeader;
	
	static std::string getPartPrefix(const std::string& msPath, size_t partIndex);
	static std::string getMetaFilename(const std::string& msPath);
};

#endif
