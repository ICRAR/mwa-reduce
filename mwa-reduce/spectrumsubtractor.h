#ifndef SPECTRUM_SUBTRACTOR_H
#define SPECTRUM_SUBTRACTOR_H

#include <boost/thread/thread.hpp>

#include <vector>

#include "lane.h"
#include "banddata.h"
#include "model.h"

namespace casa {
	class MeasurementSet;
}

class Model;
class Predicter;

class SpectrumSubtractor
{
	private:
		struct MeasureThreadInfo
		{
			casa::Complex *data;
			float *weights;
			bool *flags;
			size_t sourceIndex;
			double u, v, w;
			size_t a1, a2;
		};
		struct SubtractThreadInfo
		{
			size_t rowIndex;
			casa::Array<casa::Complex> *data;
			bool readyForWrite;
			double u, v, w;
			size_t a1, a2;
		};

	public:
		SpectrumSubtractor(casa::MeasurementSet& ms, Model& model);
		~SpectrumSubtractor();
		
		void Perform();
		void SetDataColumnName(const std::string &dataColumnName) { _dataColumnName = dataColumnName; }
		void SetFittingInterval(size_t fittingInterval) { _fittingInterval = fittingInterval; }
		const Model& RestorationModel() const { return _model; }
	private:
		void initMeasureThreadData();
		void initPredictors();
		void initSources();
		void countTimesteps(casa::ROScalarColumn<double>& timeColumn);
		
		void startMeasureThreads();
		void stopMeasureThreads();
		void measureThreadFunc(size_t threadIndex);
		
		void performSubtraction(size_t startRow, size_t endRow);
		void subtractionThreadFunc();
		void startSubtractionThreads();
		void stopSubtractionThreads();
		
		std::unique_ptr<boost::thread_group> _threadGroup;
		casa::MeasurementSet& _ms;
		Model _model;
		std::vector<ModelSource> _sources;
		BandData _bandData;
		std::vector<lane<MeasureThreadInfo>*> _taskLanes;
		size_t _cpuCount;
		std::vector<std::unique_ptr<Predicter>> _predicters;
		
		std::vector<casa::Array<casa::Complex>*> _dataBuffers;
		std::vector<casa::Array<float>*> _weightBuffers;
		std::vector<casa::Array<bool>*> _flagBuffers;
	
		std::vector<double> _spectrumSums, _spectrumWeights;

		lane<SubtractThreadInfo> _subtractWorkLane;
		lane<SubtractThreadInfo> _subtractAvailableBufferLane;
		
		static const size_t BUFFER_COUNT;
		std::string _dataColumnName;
		size_t _timestepCount, _fittingInterval;
		std::vector<double> _totalFluxPerSource, _totalFluxWeightPerSource;
		
		std::unique_ptr<casa::ArrayColumn<casa::Complex>> _dataColumn;
		std::unique_ptr<casa::ROScalarColumn<int>> _antenna1Column;
		std::unique_ptr<casa::ROScalarColumn<int>> _antenna2Column;
		std::unique_ptr<casa::ROArrayColumn<double>> _uvwColumn;
};

#endif
