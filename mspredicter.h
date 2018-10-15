#ifndef MS_PREDICTER_H
#define MS_PREDICTER_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "lane.h"
#include "predicter.h"
#include "beamevaluator.h"
#include "banddata.h"

#include "model/model.h"

#include <complex>
#include <memory>
#include <mutex>

class MSPredicter
{
public:
	struct RowData
	{
		RowData() : modelData(0)
		{	}
		
		std::complex<double> *modelData;
		size_t rowIndex, timeIndex, a1, a2;
		double u, v, w;
	};
	
	explicit MSPredicter(casacore::MeasurementSet &ms, size_t threadCount) :
		_ms(ms),
		_startChannel(0),
		_endChannel(0),
		_applyBeam(false),
		_useModelColumn(true),
		_model(),
		_laneSize(64),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
		_threadCount(threadCount)
	{ }
	
	MSPredicter(casacore::MeasurementSet &ms, size_t threadCount, const Model &model, const std::string solutionFile = "") :
		_ms(ms),
		_startChannel(0),
		_endChannel(0),
		_applyBeam(true),
		_useModelColumn(false),
		_model(model),
		_laneSize(64),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
		_threadCount(threadCount)
	{ }
	
	~MSPredicter();
	
	void Start(bool reportSources = false);
	
	bool GetNextRow(RowData& data)
	{
		return _outputLane.read(data);
	}
	void FinishRow(RowData& data)
	{
		_availableBufferLane.write(data);
	}
	
	std::mutex &IOMutex() { return _mutex; }
	
	void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
	void SetStartRow(size_t startRow) { _startRow = startRow; }
	void SetEndRow(size_t endRow) { _endRow = endRow; }
	void SetMWAPath(const std::string& mwaPath) { _mwaPath = mwaPath; }
	void SetChannelRange(size_t startChannel, size_t endChannel)
	{
		_startChannel = startChannel;
		_endChannel = endChannel;
	}

private:
	void ReadThreadFunc();
	void PredictThreadFunc();
	void clearBuffers();
	
	casacore::MeasurementSet &_ms;
	std::unique_ptr<BeamEvaluator> _beamEvaluator;
	size_t _startChannel, _endChannel;
	bool _applyBeam, _useModelColumn;
	
	Model _model;
	std::mutex _mutex;
	
	const size_t _laneSize;
	ao::lane<RowData> _workLane, _outputLane, _availableBufferLane;
	
	std::unique_ptr<std::thread> _readThread;
	std::vector<std::thread> _workThreadGroup;
	std::unique_ptr<Predicter> _predicter;
	std::vector<std::complex<double>*> _buffers;
	BandData _bandData;
	std::string _mwaPath;
	size_t _startRow, _endRow, _threadCount;
};

#endif
