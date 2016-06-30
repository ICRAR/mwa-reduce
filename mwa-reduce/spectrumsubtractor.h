#ifndef SPECTRUM_SUBTRACTOR_H
#define SPECTRUM_SUBTRACTOR_H

#include <boost/thread/thread.hpp>

#include <vector>

#include "lane.h"
#include "banddata.h"

#include "model/model.h"

namespace casacore {
	class MeasurementSet;
}

class Model;
class Predicter;

class SpectrumSubtractor
{
	private:
		struct MeasureThreadInfo
		{
			casacore::Complex *data;
			float *weights;
			bool *flags;
			size_t sourceIndex;
			double u, v, w;
			size_t a1, a2;
		};
		struct SubtractThreadInfo
		{
			size_t rowIndex;
			casacore::Array<casacore::Complex> *data;
			double u, v, w;
			size_t a1, a2;
		};

	public:
		SpectrumSubtractor(casacore::MeasurementSet& ms, Model& model);
		~SpectrumSubtractor();
		
		void Perform();
		void SetDataColumnName(const std::string &dataColumnName) { _dataColumnName = dataColumnName; }
		void SetFittingInterval(size_t fittingInterval) { _fittingInterval = fittingInterval; }
		const Model& RestorationModel() const { return _model; }
		void SetApplyBeamOnRestorationModel(bool applyBeam) { _applyBeam = applyBeam; }
	private:
		void initMeasureThreadData();
		void initPredictors();
		void initSources();
		void countTimesteps(casacore::ROScalarColumn<double>& timeColumn);
		
		void startMeasureThreads();
		void stopMeasureThreads();
		void measureThreadFunc(size_t threadIndex);
		
		void performSubtraction(size_t startRow, size_t endRow);
		void subtractionThreadFunc();
		void subtractionWriteThreadFunc();
		void startSubtractionThreads();
		void stopSubtractionThreads();
		void processWork(SubtractThreadInfo& work);
		
		std::unique_ptr<boost::thread_group> _threadGroup;
		casacore::MeasurementSet& _ms;
		Model _model;
		std::vector<ModelSource> _sources;
		BandData _bandData;
		std::vector<ao::lane<MeasureThreadInfo>*> _taskLanes;
		size_t _cpuCount;
		std::vector<std::unique_ptr<Predicter>> _predicters;
		
		std::vector<casacore::Array<casacore::Complex>*> _dataBuffers;
		std::vector<casacore::Array<float>*> _weightBuffers;
		std::vector<casacore::Array<bool>*> _flagBuffers;
	
		std::vector<double> _spectrumSums, _spectrumWeights;

		boost::mutex _subtractIOMutex;
		ao::lane<SubtractThreadInfo> _subtractWorkLane, _subtractWriteLane, _subtractAvailableLane;
		std::unique_ptr<boost::thread> _subtractWriteThread;
		
		static const size_t BUFFER_COUNT;
		std::string _dataColumnName;
		size_t _timestepCount, _fittingInterval;
		std::vector<double> _totalFluxPerSource, _totalFluxWeightPerSource;
		bool _applyBeam;
		
		std::unique_ptr<casacore::ArrayColumn<casacore::Complex>> _dataColumn;
		std::unique_ptr<casacore::ROScalarColumn<int>> _antenna1Column;
		std::unique_ptr<casacore::ROScalarColumn<int>> _antenna2Column;
		std::unique_ptr<casacore::ROArrayColumn<double>> _uvwColumn;
};

#endif
