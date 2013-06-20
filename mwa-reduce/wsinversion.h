#ifndef WS_INVERSION
#define WS_INVERSION

#include "inversionalgorithm.h"
#include "layeredimager.h"

#include "lane.h"

#include <complex>
#include <memory>

#include <casa/Arrays/Array.h>

namespace casa {
	class MeasurementSet;
}

class WSInversion : public InversionAlgorithm
{
	public:
		WSInversion() : InversionAlgorithm()
		{
		}
	
		virtual void Execute();
		
		virtual const double *ImageResult() const { return _imager->Image(); }
		virtual double ImageResultRA() const { return _phaseCentreRA; }
		virtual double ImageResultDec() const { return _phaseCentreDec; }
		virtual double ImageFrequencyHigh() const { return _freqHigh; }
		virtual double ImageFrequencyLow() const { return _freqLow; }
		virtual double ImageBeamSize() const { return _beamSize; }
	private:
		struct WorkItem
		{
			double u, v, w;
			std::complex<float> *data;
		};
		
		struct MSData
		{
			casa::MeasurementSet *ms;
			size_t channelCount, polarizationCount, matchingRows, totalRowsRead;
			double minW, maxW;
		};
		
		void initializeMeasurementSet(const std::string &measurementSet, MSData &msData);
		void gridMeasurementSet(MSData &msData);
		void countSamplesPerLayer(MSData &msData);

		void processWork(WorkItem &work)
		{
			_imager->AddData(work.data, work.u, work.v, work.w);
			delete[] work.data;
		}

		void workThread()
		{
			WorkItem workItem;
			while(_workLane->read(workItem))
			{
				processWork(workItem);
			}
		}
		void copyWeightedData(std::complex<float> *dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float> &weights, const casa::Array<bool> &flags, float rowWeight);
		void copyWeights(std::complex<float> *dest, size_t channelCount, const casa::Array<float> &weights, const casa::Array<bool> &flags, float rowWeight);
		
		std::unique_ptr<LayeredImager> _imager;
		std::unique_ptr<lane<WorkItem>> _workLane;
		double _phaseCentreRA, _phaseCentreDec;
		double _freqHigh, _freqLow;
		double _beamSize;
		double _totalWeight;
};

#endif
