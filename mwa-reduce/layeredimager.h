#ifndef LAYERED_IMAGER_H
#define LAYERED_IMAGER_H

#include "banddata.h"

#include <boost/thread/mutex.hpp>

#include <cmath>
#include <cstring>
#include <complex>
#include <vector>
#include <stack>

class LayeredImager
{
	public:
		LayeredImager(size_t width, size_t height, double pixelSizeX, double pixelSizeY, size_t fftThreadCount);
		~LayeredImager();
		
		void PrepareForObservation(size_t nWLayers, size_t maxMem, double minW, double maxW, const BandData &bandData);
		
		size_t WToLayer(double w) const
		{
			return size_t(round((fabs(w) - _minW) * (_nWLayers-1) / (_maxW - _minW)));
		}
		double LayerToW(size_t layer) const
		{
			return layer * (_maxW - _minW) / (_nWLayers-1) + _minW;
		}
		
		size_t NWLayers() const { return _nWLayers; }
		
		bool IsInLayerRange(double w) const
		{
			size_t layer = WToLayer(w);
			return layer >= layerRangeStart(_curLayerRangeIndex) && layer < layerRangeStart(_curLayerRangeIndex+1);
		}
		
		bool IsInLayerRange(double wStart, double wEnd) const
		{
			size_t
				rangeStart = layerRangeStart(_curLayerRangeIndex),
				rangeEnd = layerRangeStart(_curLayerRangeIndex+1),
				l1 = WToLayer(wStart);
			if(l1 >= rangeStart && l1 < rangeEnd)
				return true;
			size_t l2 = WToLayer(wEnd);
			return ((l2 >= rangeStart && l2 < rangeEnd) // lMax is within the range
			  || (l1 < rangeStart && l2 >= rangeEnd)  // l1 is before, l2 is after range
				|| (l2 < rangeStart && l1 >= rangeEnd) // l2 is before, l1 is after range
			);
		}
		
		size_t NPasses() const
		{
			return _nPasses;
		}
		
		void StartPass(size_t passIndex);
		
		void FinishPass();
		
		void AddData(const std::complex<float>* data, double uInM, double vInM, double wInM);
		
		void FinalizeImage(double multiplicationFactor);
		
		const double *Image() { return _imageData[0]; }
		
		size_t NFFTThreads() const { return _nFFTThreads; }
		void SetNFFTThreads(size_t nfftThreads) { _nFFTThreads = nfftThreads; }
		
	private:
		size_t layerRangeStart(size_t layerRangeIndex) const
		{
			return (_nWLayers * layerRangeIndex) / _nPasses;
		}
		void projectOnImageAndCorrect(const std::complex<double> *source, double w, size_t threadIndex);
		void initializeSqrtLMLookupTable();
		void fftThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks, size_t threadIndex);
		
		size_t _width, _height, _nWLayers, _nPasses, _curLayerRangeIndex;
		double _minW, _maxW, _pixelSizeX, _pixelSizeY;
		const BandData *_bandData;
		
		std::vector<std::vector<std::complex<double>>> _layeredUVData;
		std::vector<double*> _imageData;
		std::vector<double> _sqrtLMLookupTable;
		size_t _nFFTThreads;
};

#endif
