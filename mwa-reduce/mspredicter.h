#ifndef MS_PREDICTER_H
#define MS_PREDICTER_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include "lane.h"
#include "model.h"
#include "predicter.h"
#include "beamevaluator.h"
#include "banddata.h"

#include <boost/thread/thread.hpp>

#include <complex>
#include <memory>

class MSPredicter
{
public:
	struct RowData
	{
		RowData() : modelData(0)
		{	}
		
		std::complex<double> *modelData;
		size_t rowIndex, a1, a2;
		double u, v, w;
	};
	
	explicit MSPredicter(casa::MeasurementSet &ms, size_t threadCount) :
		_ms(ms),
		_beamEvaluator(ms, false),
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
	
	MSPredicter(casa::MeasurementSet &ms, size_t threadCount, const Model &model, const std::string solutionFile = "") :
		_ms(ms),
		_beamEvaluator(ms, false),
		_applyBeam(true),
		_useModelColumn(false),
		_model(model),
		_laneSize(64),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_solutionFile(solutionFile),
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
	
	boost::mutex &IOMutex() { return _mutex; }
	
	void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
	void SetStartRow(size_t startRow) { _startRow = startRow; }
	void SetEndRow(size_t endRow) { _endRow = endRow; }
private:
	void ReadThreadFunc();
	void PredictThreadFunc();
	void clearBuffers();
	
	casa::MeasurementSet &_ms;
	BeamEvaluator _beamEvaluator;
	size_t _channelCount;
	bool _applyBeam, _useModelColumn;
	
	Model _model;
	boost::mutex _mutex;
	
	const size_t _laneSize;
	ao::lane<RowData> _workLane, _outputLane, _availableBufferLane;
	
	std::unique_ptr<boost::thread> _readThread;
	std::unique_ptr<boost::thread_group> _workThreadGroup;
	std::unique_ptr<Predicter> _predicter;
	std::vector<std::complex<double>*> _buffers;
	std::unique_ptr<BandData> _bandData;
	std::string _solutionFile;
	size_t _startRow, _endRow, _threadCount;
};

#endif
