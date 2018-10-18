#ifndef LOFAR_MS_PREDICTER
#define LOFAR_MS_PREDICTER

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../banddata.h"
#include "../dftpredictionalgorithm.h"
#include "../lane.h"
#include "../buffered_lane.h"

#include "lbeamevaluator.h"

#include <boost/asio/io_service.hpp>
#include <boost/thread/barrier.hpp>

#include <complex>
#include <memory>
#include <thread>

class LMSPredicter
{
public:
	struct RowData
	{
		RowData() : modelData(0)
		{	}
		
		MC2x2* modelData;
		size_t rowIndex, timeIndex, a1, a2;
		double u, v, w;
	};
	
	explicit LMSPredicter(casacore::MeasurementSet &ms, size_t threadCount) :
		_ms(ms),
		_applyBeam(false),
		_useModelColumn(true),
		_dftInput(0),
		_laneSize(threadCount*16),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
		_threadCount(threadCount)
	{ }
	
	LMSPredicter(casacore::MeasurementSet &ms, size_t threadCount, DFTPredictionInput& dftInput) :
		_ms(ms),
		_applyBeam(true),
		_useModelColumn(false),
		_dftInput(&dftInput),
		_laneSize(threadCount*16),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_startRow(0),
		_endRow(ms.nrow()),
		_threadCount(threadCount)
	{ }
	
	~LMSPredicter();
	
	void Start();
	
	bool GetNextRow(RowData& data)
	{
		return _outputLane.read(data);
	}
	
	void FinishRow(RowData& data)
	{
		_availableBufferLane.write(data);
	}
	
	boost::mutex &IOMutex() { return _mutex; }
	
	void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
	void SetStartRow(size_t startRow) { _startRow = startRow; }
	void SetEndRow(size_t endRow) { _endRow = endRow; }
private:
	void ReadThreadFunc();
	void PredictThreadFunc();
	void clearBuffers();
	
	casacore::MeasurementSet &_ms;
	std::unique_ptr<LBeamEvaluator> _beamEvaluator;
	bool _applyBeam, _useModelColumn;
	
	DFTPredictionInput* _dftInput;
	boost::mutex _mutex;
	
	const size_t _laneSize;
	ao::lane<RowData> _workLane, _outputLane, _availableBufferLane;
	
	std::unique_ptr<boost::thread> _readThread;
	std::unique_ptr<boost::thread_group> _workThreadGroup;
	std::unique_ptr<DFTPredictionAlgorithm> _predicter;
	std::vector<MC2x2*> _buffers;
	BandData _bandData;
	size_t _startRow, _endRow, _threadCount;
};

#endif
