#ifndef CONTIGUOUSMS_H
#define CONTIGUOUSMS_H

#include "msprovider.h"

#include "../msselection.h"
#include "../multibanddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <memory>

class ContiguousMS : public MSProvider
{
public:
	ContiguousMS(const string& msPath, const std::string& dataColumnName, MSSelection selection, PolarizationEnum polOut, bool includeModel);
	
	virtual casa::MeasurementSet &MS() { return _ms; }
	
	virtual size_t RowId() const { return _row; }
	
	virtual bool NextRow();
	
	virtual void Reset();
	
	virtual void ReadMeta(double& u, double& v, double& w, size_t& dataDescId);
	
	virtual void ReadData(std::complex<float>* buffer);
	
	virtual void ReadModel(std::complex<float>* buffer);
	
	virtual void WriteModel(size_t rowId, std::complex<float>* buffer);
	
	virtual void ReadWeights(float* buffer);
	
	virtual void ReadWeights(std::complex<float>* buffer);
	
	virtual void ReopenRW()
	{
		_ms.reopenRW();
	}
private:
	size_t _row;
	size_t _timestep;
	double _time;
	int _dataDescId;
	bool _isMetaRead, _isDataRead, _isModelRead, _isWeightRead;
	bool _isModelColumnPrepared;
	size_t _startRow, _endRow;
	std::vector<PolarizationEnum> _inputPolarizations;
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
	
	void prepareModelColumn();
	void readMeta()
	{
		if(!_isMetaRead)
		{
			_dataDescId = _dataDescIdColumn(_row);
			_isMetaRead = true;
		}
	}
	void readData()
	{
		if(!_isDataRead)
		{
			_dataColumn.get(_row, _dataArray);
			_isDataRead = true;
		}
	}
	void readWeights()
	{
		if(!_isWeightRead)
		{
			_flagColumn.get(_row, _flagArray);
			if(_msHasWeights)
				_weightColumn.get(_row, _weightArray);
			_isWeightRead = true;
		}
	}
	void readModel()
	{
		if(!_isModelRead)
		{
			_modelColumn->get(_row, _modelArray);
			_isModelRead = true;
		}
	}
};

#endif
