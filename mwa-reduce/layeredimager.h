#ifndef LAYERED_IMAGER_H
#define LAYERED_IMAGER_H

#include "banddata.h"
#include "multibanddata.h"

#include <boost/thread/mutex.hpp>

#include <cmath>
#include <cstring>
#include <complex>
#include <vector>
#include <stack>

template<typename NumType>
class ImageBufferAllocator;

class LayeredImager
{
	public:
		enum GridModeEnum { NearestNeighbour, KaiserBessel };
		
		typedef double imgnum_t;
		
		LayeredImager(size_t width, size_t height, double pixelSizeX, double pixelSizeY, size_t fftThreadCount, ImageBufferAllocator<double>* allocator, size_t kernelSize = 7, size_t overSamplingFactor = 63);
		~LayeredImager();
		
		void PrepareWLayers(size_t nWLayers, double maxMem, double minW, double maxW);
		
		void PrepareBand(const MultiBandData &bandData)
		{
			_bandData = bandData;
		}
		
		size_t WToLayer(double wInLambda) const
		{
			if(_nWLayers == 1)
				return 0;
			else {
				if(_isComplex)
					return size_t(round((wInLambda + _maxW) * (_nWLayers-1) / (_maxW + _maxW)));
				else
					return size_t(round((fabs(wInLambda) - _minW) * (_nWLayers-1) / (_maxW - _minW)));
			}
		}
		
		double LayerToW(size_t layer) const
		{
			if(_nWLayers == 1)
				return 0.0;
			else {
				if(_isComplex)
					return layer * (_maxW + _maxW) / (_nWLayers-1) - _maxW;
				else
					return layer * (_maxW - _minW) / (_nWLayers-1) + _minW;
			}
		}
		
		size_t NWLayers() const { return _nWLayers; }
		
		bool IsInLayerRange(double wInLambda) const
		{
			size_t layer = WToLayer(wInLambda);
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
		
		// Inversion (uv to image) methods
		void AddData(const std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM);
		
		void AddDataSample(std::complex<float> sample, double uInLambda, double vInLambda, double wInLambda);
		
		void StartInversionPass(size_t passIndex);
		
		void FinishInversionPass();
		
		void FinalizeImage(double multiplicationFactor);
		
		// Sampling (image to uv) methods
		void InitializePrediction(const double *image, double multiplicationFactor)
		{
			initializePrediction(image, multiplicationFactor, _imageData);
		}
		
		void InitializePrediction(const double *real, const double *imaginary, double multiplicationFactor)
		{
			initializePrediction(real, multiplicationFactor, _imageData);
			initializePrediction(imaginary, multiplicationFactor, _imageDataImaginary);
		}

		void StartPredictionPass(size_t passIndex);
		
		void SampleData(std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM);
		
		double *RealImage() { return _imageData[0]; }
		double *ImaginaryImage() { return _imageDataImaginary[0]; }
		
		size_t NFFTThreads() const { return _nFFTThreads; }
		void SetNFFTThreads(size_t nfftThreads) { _nFFTThreads = nfftThreads; }
		
		enum GridModeEnum GridMode() const { return _gridMode; }
		void SetGridMode(enum GridModeEnum mode) { _gridMode = mode; }
		
		void SetIsComplex(bool isComplex) { _isComplex = isComplex; }
		
		void SetImageConjugatePart(bool imageConjugatePart) { _imageConjugatePart = imageConjugatePart; }
		
		void SetDenormalPhaseCentre(double dl, double dm) { _phaseCentreDL = dl; _phaseCentreDM = dm; }
		
		void GetGriddingCorrectionImage(double* image) const;
		
		void ReplaceRealImageBuffer(double* newBuffer);
		
		void ReplaceImaginaryImageBuffer(double* newBuffer);
	private:
		size_t layerRangeStart(size_t layerRangeIndex) const
		{
			return (_nWLayers * layerRangeIndex) / _nPasses;
		}
		template<bool IsComplex>
		void projectOnImageAndCorrect(const std::complex<double> *source, double w, size_t threadIndex);
		template<bool IsComplex>
		void copyImageToLayerAndInverseCorrect(std::complex<double> *dest, double w);
		void initializeSqrtLMLookupTable();
		void initializeSqrtLMLookupTableForSampling();
		void initializeLayeredUVData(size_t n);
		void freeLayeredUVData() { initializeLayeredUVData(0); }
		void fftToImageThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks, size_t threadIndex);
		void fftToUVThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks);
		void finalizeImage(double multiplicationFactor, std::vector<double*>& dataArray);
		void initializePrediction(const double *image, double multiplicationFactor, std::vector<double*>& dataArray);
		
		void makeKernels();
		void makeKernel(std::vector<double> &kernel, double alpha, size_t overSamplingFactor);
		double bessel0(double x, double precision);
		template<bool Inverse>
		void correctImageForKernel(double *image) const;
		
		size_t _width, _height, _nWLayers, _nPasses, _curLayerRangeIndex;
		double _minW, _maxW, _pixelSizeX, _pixelSizeY, _phaseCentreDL, _phaseCentreDM;
		bool _isComplex, _imageConjugatePart;
		MultiBandData _bandData;
		
		enum GridModeEnum _gridMode;
		size_t _overSamplingFactor, _kernelSize;
		std::vector<double> _1dKernel;
		std::vector<std::vector<double>> _griddingKernels;
		
		std::vector<std::complex<double>*> _layeredUVData;
		std::vector<double*> _imageData, _imageDataImaginary;
		std::vector<double> _sqrtLMLookupTable;
		size_t _nFFTThreads;
		ImageBufferAllocator<double>* _imageBufferAllocator;
};

#endif
