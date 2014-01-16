#ifndef WS_INVERSION
#define WS_INVERSION

#include "inversionalgorithm.h"
#include "layeredimager.h"

#include "lane.h"
#include "multibanddata.h"

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
			size_t dataDescId;
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
			size_t rowIndex, dataDescId;
		};
		
		struct MSData
		{
			public:
				MSData();
				~MSData();
				std::unique_ptr<casa::MeasurementSet> ms;
				MultiBandData bandData;
				size_t startChannel, endChannel;
				size_t polarizationCount, matchingRows, totalRowsProcessed;
				double minW, maxW;
				size_t rowStart, rowEnd;
			
				MultiBandData SelectedBand() const { return MultiBandData(bandData, startChannel, endChannel); }
			private:
				MSData(const MSData &source);
				
				void operator=(const MSData &source);
		};
		
		void initializeMeasurementSet(const std::string &measurementSet, MSData &msData);
		void gridMeasurementSet(MSData &msData);
		void countSamplesPerLayer(MSData &msData);

		void sampleToMeasurementSet(MSData &msData);

		void workThread(ao::lane<InversionWorkItem>* workLane)
		{
			InversionWorkItem workItem;
			while(workLane->read(workItem))
			{
				_imager->AddData(workItem.data, workItem.dataDescId, workItem.u, workItem.v, workItem.w);
				delete[] workItem.data;
			}
		}
		
		void workThreadParallel(const MultiBandData* selectedBand);
		void workThreadPerSample(ao::lane<InversionWorkSample>* workLane);
		
		void visSampleCalcThread(ao::lane<SamplingWorkItem>* inputLane, ao::lane<SamplingWorkItem>* outputLane);
		void visSampleWriteThread(ao::lane<SamplingWorkItem>* samplingWorkLane, const MSData* msData);
		void copyWeightedData(std::complex<float> *dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float> &weights, const casa::Array<bool> &flags, float rowWeight);
		void copyWeights(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight);
		int polarizationIndex(int polCount) const
		{
			if(polCount == 4)
			{
				switch(Polarization())
				{
					default: return 0;
					case Polarization::XY: return 1;
					case Polarization::YX: return 2;
					case Polarization::YY: return 3;
				}
			}
			else if(polCount == 2)
			{
				switch(Polarization())
				{
					default: return 0;
					case Polarization::YY: return 1;
					case Polarization::XY:
					case Polarization::YX:
						throw std::runtime_error("Invalid polarization");
				}
			}
			else
				return 0;
		}

		std::unique_ptr<LayeredImager> _imager;
		std::unique_ptr<ao::lane<InversionWorkItem>> _inversionWorkLane;
		std::unique_ptr<casa::ArrayColumn<casa::Complex>> _modelColumn;
		double _phaseCentreRA, _phaseCentreDec;
		bool _hasFrequencies;
		double _freqHigh, _freqLow;
		double _bandStart, _bandEnd;
		double _beamSize;
		double _totalWeight;
		double _startTime;
		LayeredImager::GridModeEnum _gridMode;
		size_t _cpuCount, _laneBufferSize;
		int64_t _memSize;
};

#endif
