#ifndef WS_INVERSION
#define WS_INVERSION

#include "inversionalgorithm.h"
#include "layeredimager.h"

#include "lane.h"

#include <complex>
#include <memory>

#include <casa/Arrays/Array.h>
#include <tables/Tables/ArrayColumn.h>

namespace casa {
	class MeasurementSet;
}

class WSInversion : public InversionAlgorithm
{
	public:
		WSInversion();
	
		virtual void Invert();
		
		virtual void InvertToVisibilities(const double *image);
		
		virtual const double *ImageResult() const { return _imager->Image(); }
		virtual double ImageResultRA() const { return _phaseCentreRA; }
		virtual double ImageResultDec() const { return _phaseCentreDec; }
		virtual double ImageHighestFrequencyChannel() const { return _freqHigh; }
		virtual double ImageLowestFrequencyChannel() const { return _freqLow; }
		virtual double ImageBandStart() const { return _bandStart; }
		virtual double ImageBandEnd() const { return _bandEnd; }
		virtual double ImageBeamSize() const { return _beamSize; }
		virtual double ImageStartTime() const { return _startTime; }
		
		enum LayeredImager::GridModeEnum GridMode() const { return _gridMode; }
		void SetGridMode(LayeredImager::GridModeEnum gridMode) { _gridMode = gridMode; }
		
		virtual bool HasGriddingCorrectionImage() const { return _gridMode == LayeredImager::KaiserBessel; }
		virtual void GetGriddingCorrectionImage(double *image) const { _imager->GetGriddingCorrectionImage(image); }
	private:
		struct InversionWorkItem
		{
			double u, v, w;
			std::complex<float> *data;
		};
		struct InversionWorkSample
		{
			double uInLambda, vInLambda, wInLambda;
			std::complex<float> sample;
		};
		struct SamplingWorkItem
		{
			double u, v, w;
			std::complex<float> *data;
			size_t rowIndex;
		};
		
		struct MSData
		{
			public:
				MSData();
				~MSData();
				std::unique_ptr<casa::MeasurementSet> ms;
				size_t channelCount, polarizationCount, matchingRows, totalRowsProcessed;
				double minW, maxW;
				size_t rowStart, rowEnd;
			
			private:
				MSData(const MSData &source);
				
				void operator=(const MSData &source);
		};
		
		void initializeMeasurementSet(const std::string &measurementSet, MSData &msData);
		void gridMeasurementSet(MSData &msData);
		void countSamplesPerLayer(MSData &msData);

		void sampleToMeasurementSet(MSData &msData);

		void workThread(lane<InversionWorkItem>* workLane)
		{
			InversionWorkItem workItem;
			while(workLane->read(workItem))
			{
				_imager->AddData(workItem.data, workItem.u, workItem.v, workItem.w);
				delete[] workItem.data;
			}
		}
		
		void workThreadParallel(BandData* bandData);
		void workThreadPerSample(lane<InversionWorkSample>* workLane);
		
		void visSampleCalcThread(lane<SamplingWorkItem>* inputLane, lane<SamplingWorkItem>* outputLane);
		void visSampleWriteThread(lane<SamplingWorkItem>* samplingWorkLane);
		void copyWeightedData(std::complex<float> *dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float> &weights, const casa::Array<bool> &flags, float rowWeight);
		void copyWeights(std::complex<float>* dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight);
		int polarizationIndex() const
		{
			switch(Polarization())
			{
				default: return 0;
				case XY: return 1;
				case YX: return 2;
				case YY: return 3;
			}
		}

		std::unique_ptr<LayeredImager> _imager;
		std::unique_ptr<lane<InversionWorkItem>> _inversionWorkLane;
		std::unique_ptr<casa::ArrayColumn<casa::Complex>> _modelColumn;
		double _phaseCentreRA, _phaseCentreDec;
		bool _hasFrequencies;
		double _freqHigh, _freqLow;
		double _bandStart, _bandEnd;
		double _beamSize;
		double _totalWeight;
		double _startTime;
		LayeredImager::GridModeEnum _gridMode;
		size_t _cpuCount;
};

#endif
