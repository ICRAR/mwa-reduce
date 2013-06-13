#ifndef WS_INVERSION
#define WS_INVERSION

#include "inversionalgorithm.h"
#include "layeredimager.h"

#include "lane.h"

#include <complex>
#include <memory>

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
	private:
		struct WorkItem
		{
			double u, v, w;
			std::complex<float> *data;
		};

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
		
		std::unique_ptr<LayeredImager> _imager;
		std::unique_ptr<lane<WorkItem>> _workLane;
		double _phaseCentreRA, _phaseCentreDec;
};

#endif
