#include "layeredimager.h"
#include "imagebufferallocator.h"

#include <fftw3.h>

#include <iostream>
#include <fstream>
#include <stack>

#include <boost/thread/thread.hpp>

LayeredImager::LayeredImager(size_t width, size_t height, double pixelSizeX, double pixelSizeY, size_t fftThreadCount, ImageBufferAllocator<double>* allocator, size_t kernelSize, size_t overSamplingFactor) :
	_width(width),
	_height(height),
	_pixelSizeX(pixelSizeX),
	_pixelSizeY(pixelSizeY),
	_phaseCentreDL(0.0),
	_phaseCentreDM(0.0),
	_isComplex(false),
	_imageConjugatePart(false),
	_gridMode(NearestNeighbour),
	_overSamplingFactor(overSamplingFactor),
	_kernelSize(kernelSize),
	_imageData(fftThreadCount),
	_imageDataImaginary(fftThreadCount),
	_nFFTThreads(fftThreadCount),
	_imageBufferAllocator(allocator)
{
	makeKernels();
}

LayeredImager::~LayeredImager()
{
	for(size_t i=0; i!=_nFFTThreads; ++i)
	{
		_imageBufferAllocator->Free(_imageData[i]);
		_imageBufferAllocator->Free(_imageDataImaginary[i]);
	}
	freeLayeredUVData();
	fftw_cleanup();
}

void LayeredImager::PrepareWLayers(size_t nWLayers, double maxMem, double minW, double maxW)
{
	_minW = minW;
	_maxW = maxW;
	_nWLayers = nWLayers;
	
	size_t nrCopies = _nFFTThreads;
	if(nrCopies > _nWLayers) nrCopies = _nWLayers;
	double memPerImage = _width * _height * sizeof(double);
	double memPerCore = memPerImage * 5.0; // two complex ones for FFT, one for projecting on
	double remainingMem = maxMem - nrCopies * memPerCore;
	if(remainingMem <= memPerImage * _nFFTThreads)
	{
		_nFFTThreads = size_t(maxMem*3.0/(5.0*memPerCore)); // times 3/5 to use 3/5 of mem for FFTing at most
		if(_nFFTThreads==0) _nFFTThreads = 1;
		remainingMem = maxMem - _nFFTThreads * memPerCore;
		
		std::cout <<
			"WARNING: the amount of available memory is too low for the image size,\n"
			"       : not all cores might be used.\n"
			"       : nr buffers avail for FFT: " << _nFFTThreads << " remaining mem: " << round(remainingMem/1.0e8)/10.0 << " GB \n";
	}
	
	// Allocate FFT buffers
	size_t imgSize = _height * _width;
	for(size_t i=0; i!=_nFFTThreads; ++i)
	{
		_imageData[i] = _imageBufferAllocator->Allocate(imgSize);
		memset(_imageData[i], 0, imgSize * sizeof(double));
		if(_isComplex)
		{
			_imageDataImaginary[i] = _imageBufferAllocator->Allocate(imgSize);
			memset(_imageDataImaginary[i], 0, imgSize * sizeof(double));
		}
	}
	
	// Calculate nr wlayers per pass from remaining memory
	int maxNWLayersPerPass = int((double) remainingMem / (2.0*memPerImage));
	if(maxNWLayersPerPass < 1)
		maxNWLayersPerPass=1;
	_nPasses = (nWLayers+maxNWLayersPerPass-1)/maxNWLayersPerPass;
	if(_nPasses == 0) _nPasses = 1;
	std::cout << "Will process " << (_nWLayers / _nPasses) << "/" << _nWLayers << " w-layers per pass.\n";
	
	_curLayerRangeIndex = 0;
}

void LayeredImager::initializeLayeredUVData(size_t n)
{
	while(_layeredUVData.size() > n)
	{
		_imageBufferAllocator->Free(_layeredUVData.back());
		_layeredUVData.pop_back();
	}
	while(_layeredUVData.size() < n)
		_layeredUVData.push_back(_imageBufferAllocator->AllocateComplex(_width * _height));
}

void LayeredImager::StartInversionPass(size_t passIndex)
{
	initializeSqrtLMLookupTable();
	
	_curLayerRangeIndex = passIndex;
	size_t nLayersInPass = layerRangeStart(passIndex+1) - layerRangeStart(passIndex);
	initializeLayeredUVData(nLayersInPass);
	for(size_t i=0; i!=nLayersInPass; ++i)
		memset(_layeredUVData[i], 0, _width*_height * sizeof(double)*2);
}

void LayeredImager::StartPredictionPass(size_t passIndex)
{
	initializeSqrtLMLookupTableForSampling();
	
	_curLayerRangeIndex = passIndex;
	size_t layerOffset = layerRangeStart(passIndex);
	size_t nLayersInPass = layerRangeStart(passIndex+1) - layerOffset;
	initializeLayeredUVData(nLayersInPass);
	
	std::stack<size_t> layers;
	for(size_t layer=0; layer!=nLayersInPass; ++layer)
		layers.push(layer);
	
	boost::mutex mutex;
	boost::thread_group threadGroup;
	for(size_t i=0; i!=_nFFTThreads; ++i)
		threadGroup.add_thread(new boost::thread(&LayeredImager::fftToUVThreadFunction, this, &mutex, &layers));
	threadGroup.join_all();
}

void LayeredImager::fftToImageThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks, size_t threadIndex)
{
	const size_t imgSize = _width * _height;
	std::complex<double> *fftwIn = _imageBufferAllocator->AllocateComplex(imgSize);
	// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	std::complex<double> *fftwOut = _imageBufferAllocator->AllocateComplex(imgSize);
	// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	
	boost::mutex::scoped_lock lock(*mutex);
	fftw_plan plan =
		fftw_plan_dft_2d(_width, _height,
			reinterpret_cast<fftw_complex*>(fftwIn), reinterpret_cast<fftw_complex*>(fftwOut),
			FFTW_BACKWARD, FFTW_ESTIMATE);
		
	const size_t layerOffset = layerRangeStart(_curLayerRangeIndex);

	while(!tasks->empty())
	{
		size_t layer = tasks->top();
		tasks->pop();
		lock.unlock();
		
		// Fourier transform the layer
		std::complex<double> *uvData = _layeredUVData[layer];
		memcpy(fftwIn, uvData, imgSize * sizeof(double) * 2);
		fftw_execute(plan);
		
		// Add layer to full image
		if(_isComplex)
			projectOnImageAndCorrect<true>(fftwOut, LayerToW(layer + layerOffset), threadIndex);
		else
			projectOnImageAndCorrect<false>(fftwOut, LayerToW(layer + layerOffset), threadIndex);
		
		// lock for accessing tasks in guard
		lock.lock();
	}
	lock.unlock();
	fftw_destroy_plan(plan);
	_imageBufferAllocator->Free(fftwIn);
	_imageBufferAllocator->Free(fftwOut);
}

void LayeredImager::fftToUVThreadFunction(boost::mutex *mutex, std::stack<size_t> *tasks)
{
	const size_t imgSize = _width * _height;
	std::complex<double> *fftwIn = _imageBufferAllocator->AllocateComplex(imgSize);
		// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	std::complex<double> *fftwOut = _imageBufferAllocator->AllocateComplex(imgSize);
		// reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	
	boost::mutex::scoped_lock lock(*mutex);
	fftw_plan plan =
		fftw_plan_dft_2d(_width, _height,
			reinterpret_cast<fftw_complex*>(fftwIn), reinterpret_cast<fftw_complex*>(fftwOut),
			FFTW_FORWARD, FFTW_ESTIMATE);
		
	const size_t layerOffset = layerRangeStart(_curLayerRangeIndex);

	while(!tasks->empty())
	{
		size_t layer = tasks->top();
		tasks->pop();
		lock.unlock();
		
		// Make copy of input and w-correct it
		if(_isComplex)
			copyImageToLayerAndInverseCorrect<true>(fftwIn, LayerToW(layer + layerOffset));
		else
			copyImageToLayerAndInverseCorrect<false>(fftwIn, LayerToW(layer + layerOffset));
		
		// Fourier transform the layer
		fftw_execute(plan);
		std::complex<double> *uvData = _layeredUVData[layer];
		memcpy(uvData, fftwOut, imgSize * sizeof(double) * 2);
		
		// lock for accessing tasks in guard
		lock.lock();
	}
	lock.unlock();
	
	_imageBufferAllocator->Free(fftwIn);
	_imageBufferAllocator->Free(fftwOut);
}

void LayeredImager::FinishInversionPass()
{
	size_t layerOffset = layerRangeStart(_curLayerRangeIndex);
	size_t nPlanes = layerRangeStart(_curLayerRangeIndex+1) - layerOffset;
	std::stack<size_t> planes;
	for(size_t plane=0; plane!=nPlanes; ++plane)
		planes.push(plane);
	
	boost::mutex mutex;
	boost::thread_group threadGroup;
	for(size_t i=0; i!=_nFFTThreads; ++i)
		threadGroup.add_thread(new boost::thread(&LayeredImager::fftToImageThreadFunction, this, &mutex, &planes, i));
	threadGroup.join_all();
}

void LayeredImager::makeKernels()
{
	_griddingKernels.resize(_overSamplingFactor * _overSamplingFactor);
	_1dKernel.resize(_kernelSize*_overSamplingFactor);
	const double alpha = _kernelSize;
	makeKernel(_1dKernel, alpha, _overSamplingFactor);
	
	std::vector<std::vector<double>>::iterator gridKernelIter = _griddingKernels.begin();
	for(size_t j=0; j!=_overSamplingFactor; ++j)
	{
		for(size_t i=0; i!=_overSamplingFactor; ++i)
		{
			std::vector<double> &kernel = _griddingKernels[(_overSamplingFactor-j-1)*_overSamplingFactor + _overSamplingFactor - i - 1];
			kernel.resize(_kernelSize * _kernelSize);
			std::vector<double>::iterator kernelValueIter = kernel.begin();
			for(size_t y=0; y!=_kernelSize; ++y)
			{
				for(size_t x=0; x!=_kernelSize; ++x)
				{
					size_t xIndex = x*_overSamplingFactor + i, yIndex = y*_overSamplingFactor + j;
					*kernelValueIter = _1dKernel[xIndex] * _1dKernel[yIndex];
					++kernelValueIter;
				}
			}
			++gridKernelIter;
		}
	}
}

void LayeredImager::makeKernel(std::vector<double> &kernel, double alpha, size_t overSamplingFactor)
{
	size_t
		n = kernel.size(),
		mid = n/2;
	std::vector<double> sincKernel(mid+1);
	const double filterRatio = 1.0 / double(overSamplingFactor); // FILTER POINT / TOTAL BANDWIDTH
	sincKernel[0] = filterRatio;
	for(size_t i=1; i!=mid+1; i++)
	{
		double x = i;
		sincKernel[i] = sin(M_PI*filterRatio*x)/(M_PI*x);
	}
	const double
		normFactor = double(overSamplingFactor) / bessel0(alpha, 1e-8),
		nSquared = n*n;
	for(size_t i=0; i!=mid+1; i++)
		kernel[mid+i] = sincKernel[i] * bessel0(alpha * sqrt(1.0-(double(i*i)/nSquared)), 1e-8) * normFactor;
	for(size_t i=0; i!=mid; i++)
		kernel[i] = kernel[n-1-i];
}

double LayeredImager::bessel0(double x, double precision)
{
	// Calculate I_0 = SUM of m 0 -> inf [ (x/2)^(2m) ]
	// This is the unnormalized bessel function of order 0.
	double
		d = 0.0,
		ds = 1.0,
		s = 1.0;
	do
	{
		d += 2.0;
		ds *= x*x/(d*d);
		s += ds;
	} while (ds > s*precision);
	return s;
}

void LayeredImager::AddData(const std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM)
{
	const BandData& curBand = _bandData[dataDescId];
	for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
	{
		double
			wavelength = curBand.ChannelWavelength(ch),
			u = uInM / wavelength,
			v = vInM / wavelength,
			w = wInM / wavelength;
		AddDataSample(data[ch], u, v, w);
	}
}

void LayeredImager::AddDataSample(std::complex<float> sample, double uInLambda, double vInLambda, double wInLambda)
{
 	const size_t
		layerOffset = layerRangeStart(_curLayerRangeIndex),
		layerRangeEnd = layerRangeStart(_curLayerRangeIndex+1);
	if(_imageConjugatePart)
	{
		uInLambda = -uInLambda;
		vInLambda = -vInLambda;
		sample = std::conj(sample);
	}
	if(wInLambda < 0.0 && !_isComplex)
	{
		uInLambda = -uInLambda;
		vInLambda = -vInLambda;
		wInLambda = -wInLambda;
		sample = std::conj(sample);
	}
	size_t
		wLayer = WToLayer(wInLambda);
	if(wLayer >= layerOffset && wLayer < layerRangeEnd)
	{
		size_t layerIndex = wLayer - layerOffset;
		std::complex<double>* uvData = _layeredUVData[layerIndex];
		if(_gridMode == KaiserBessel)
		{
			double
				xExact = uInLambda * _pixelSizeX * _width,
				yExact = vInLambda * _pixelSizeY * _height;
			int
				x = round(xExact),
				y = round(yExact),
				xKernel = round((xExact - double(x)) * _overSamplingFactor),
				yKernel = round((yExact - double(y)) * _overSamplingFactor);
			xKernel = (xKernel + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			yKernel = (yKernel + (_overSamplingFactor*3)/2) % _overSamplingFactor;
			std::vector<double> &kernel = _griddingKernels[xKernel + yKernel*_overSamplingFactor];
			int mid = _kernelSize / 2;
			if(x < 0) x += _width;
			if(y < 0) y += _height;
			if(x >= 0 && y >= 0 && x < (int) _width && y < (int) _height)
			{
				// Are we on the edge?
				if(x < mid || x+mid+1 >= int(_width) || y < mid || y+mid+1 >= int(_height))
				{
					std::vector<double>::iterator kernelIter = kernel.begin();
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						size_t cy = ((y+j+_height-mid) % _height) * _width;
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							size_t cx = (x+i+_width-mid) % _width;
							std::complex<double> *uvRowPtr = &uvData[cx + cy];
							*uvRowPtr += std::complex<double>(sample.real() * (*kernelIter), sample.imag() * (*kernelIter));
							++kernelIter;
						}
					}
				}
				else {
					x -= mid;
					y -= mid;
					std::vector<double>::iterator kernelIter = kernel.begin();
					for(size_t j=0; j!=_kernelSize; ++j)
					{
						std::complex<double> *uvRowPtr = &uvData[x + y*_width];
						for(size_t i=0; i!=_kernelSize; ++i)
						{
							*uvRowPtr += std::complex<double>(sample.real() * (*kernelIter), sample.imag() * (*kernelIter));
							++uvRowPtr;
							++kernelIter;
						}
						++y;
					}
				}
			}
		}
		else {
			int
				x = int(round(uInLambda * _pixelSizeX * _width)),
				y = int(round(vInLambda * _pixelSizeY * _height));
			if(x < 0) x += _width;
			if(y < 0) y += _height;
			if(x >= 0 && y >= 0 && x < (int) _width && y < (int) _height)
			{
				uvData[x + y*_width] += sample;
			} else {
				//std::cout << "Sample fell off uv-plane (" << x << "," << y << ")\n";
			}
		}
	}
}

void LayeredImager::SampleData(std::complex<float>* data, size_t dataDescId, double uInM, double vInM, double wInM)
{
 	const size_t
		layerOffset = layerRangeStart(_curLayerRangeIndex),
		layerRangeEnd = layerRangeStart(_curLayerRangeIndex+1);
	const BandData& curBand(_bandData[dataDescId]);
	for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
	{
		double
			wavelength = curBand.ChannelWavelength(ch),
			u = uInM / wavelength,
			v = vInM / wavelength,
			w = wInM / wavelength;
		bool isConjugated = (w < 0.0);
		if(isConjugated)
		{
			u = -u;
			v = -v;
			w = -w;
		}
		size_t
			wLayer = WToLayer(w);
		if(wLayer >= layerOffset && wLayer < layerRangeEnd)
		{
			size_t layerIndex = wLayer - layerOffset;
			std::complex<double>* uvData = _layeredUVData[layerIndex];
			int
				x = int(round(u * _pixelSizeX * _width)),
				y = int(round(v * _pixelSizeY * _height));
			if(x < 0) x += _width;
			if(y < 0) y += _height;
			std::complex<double> sample;
			if(_gridMode == KaiserBessel)
			{
				sample = 0.0;
				double
					xExact = u * _pixelSizeX * _width,
					yExact = v * _pixelSizeY * _height;
				int
					x = round(xExact),
					y = round(yExact),
					xKernel = round((xExact - double(x)) * _overSamplingFactor),
					yKernel = round((yExact - double(y)) * _overSamplingFactor);
				xKernel = (xKernel + (_overSamplingFactor*3)/2) % _overSamplingFactor;
				yKernel = (yKernel + (_overSamplingFactor*3)/2) % _overSamplingFactor;
				std::vector<double> &kernel = _griddingKernels[xKernel + yKernel*_overSamplingFactor];
				int mid = _kernelSize / 2;
				if(x < 0) x += _width;
				if(y < 0) y += _height;
				if(x >= 0 && y >= 0 && x < (int) _width && y < (int) _height)
				{
					// Are we on the edge?
					if(x < mid || x+mid+1 >= int(_width) || y < mid || y+mid+1 >= int(_height))
					{
						std::vector<double>::iterator kernelIter = kernel.begin();
						for(size_t j=0; j!=_kernelSize; ++j)
						{
							size_t cy = ((y+j+_height-mid) % _height) * _width;
							for(size_t i=0; i!=_kernelSize; ++i)
							{
								size_t cx = (x+i+_width-mid) % _width;
								std::complex<double> *uvRowPtr = &uvData[cx + cy];
								sample += std::complex<double>(uvRowPtr->real() * (*kernelIter), uvRowPtr->imag() * (*kernelIter));
								++kernelIter;
							}
						}
					}
					else {
						x -= mid;
						y -= mid;
						std::vector<double>::iterator kernelIter = kernel.begin();
						for(size_t j=0; j!=_kernelSize; ++j)
						{
							std::complex<double> *uvRowPtr = &uvData[x + y*_width];
							for(size_t i=0; i!=_kernelSize; ++i)
							{
								sample += std::complex<double>(uvRowPtr->real() * (*kernelIter), uvRowPtr->imag() * (*kernelIter));
								++uvRowPtr;
								++kernelIter;
							}
							++y;
						}
					}
				}
			}
			else {
				if(x >= 0 && y >= 0 && x < (int) _width && y < (int) _height)
				{
					sample = uvData[x + y*_width];
				} else {
					//std::cout << "Sampling outside uv-plane (" << x << "," << y << ")\n";
				}
			}
			if(isConjugated)
				data[ch] = sample;
			else
				data[ch] = std::conj(sample);
		} else {
			data[ch] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
		}
	}
}

void LayeredImager::FinalizeImage(double multiplicationFactor)
{
	freeLayeredUVData();
	finalizeImage(multiplicationFactor, _imageData);
	if(_isComplex)
		finalizeImage(multiplicationFactor, _imageDataImaginary);
}

void LayeredImager::finalizeImage(double multiplicationFactor, std::vector<double*>& dataArray)
{
	for(size_t i=1;i!=_nFFTThreads;++i)
	{
		double *primaryData = dataArray[0];
		double *endPtr = dataArray[i] + (_width * _height);
		for(double *dataPtr = dataArray[i]; dataPtr!=endPtr; ++dataPtr)
		{
			*primaryData += *dataPtr;
			++primaryData;
		}
	}
	double *dataPtr = dataArray[0];
	for(size_t y=0;y!=_height;++y)
	{
		double m = ((double) y-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		for(size_t x=0;x!=_width;++x)
		{
			double l = ((_width/2)-(double) x) * _pixelSizeX + _phaseCentreDL;
			*dataPtr *= multiplicationFactor * sqrt(1.0 - l*l - m*m);
			++dataPtr;
		}
	}
	
	if(_gridMode == KaiserBessel)
		correctImageForKernel<false>(dataArray[0]);
}

template<bool Inverse>
void LayeredImager::correctImageForKernel(double *image) const
{
	const size_t n = _width * _overSamplingFactor;
	double
		*fftwIn = reinterpret_cast<double*>(fftw_malloc(n/2 * sizeof(double))),
		*fftwOut = reinterpret_cast<double*>(fftw_malloc(n/2 * sizeof(double)));
	
	fftw_plan plan = fftw_plan_r2r_1d(n/2, fftwIn, fftwOut, FFTW_REDFT01, FFTW_ESTIMATE);
	memset(fftwIn, 0, n/2 * sizeof(double));
	memcpy(fftwIn, &_1dKernel[_kernelSize*_overSamplingFactor/2], (_kernelSize*_overSamplingFactor/2+1) * sizeof(double));
	fftw_execute(plan);
	fftw_free(fftwIn);
	
	double normFactor = 1.0 / (_overSamplingFactor * _overSamplingFactor);
	for(size_t y=0; y!=_height; ++y)
	{
		for(size_t x=0; x!=_width; ++x)
		{
			double xVal = (x>=_width/2) ? fftwOut[x-_width/2] : fftwOut[_width/2-x];
			double yVal = (y>=_height/2) ? fftwOut[y-_height/2] : fftwOut[_height/2-y];
			if(Inverse)
				*image *= xVal * yVal * normFactor;
			else
				*image /= xVal * yVal * normFactor;
			++image;
		}
	}
	
	fftw_destroy_plan(plan);
	fftw_free(fftwOut);
}

void LayeredImager::GetGriddingCorrectionImage(double *image) const
{
	for(size_t i=0; i!=_width*_height; ++i)
		image[i] = 1.0;
	correctImageForKernel<true>(image);
}

void LayeredImager::initializePrediction(const double* image, double multiplicationFactor, std::vector<double*>& dataArray)
{
	double *dataPtr = dataArray[0];
	const double *inPtr = image;
	for(size_t y=0;y!=_height;++y)
	{
		double m = ((double) y-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		for(size_t x=0;x!=_width;++x)
		{
			double l = ((_width/2)-(double) x) * _pixelSizeX + _phaseCentreDL;
			if(std::isfinite(*dataPtr) && l*l + m*m < 1.0)
				*dataPtr = *inPtr * multiplicationFactor / sqrt(1.0 - l*l - m*m);
			else
				*dataPtr = 0.0;
			++dataPtr;
			++inPtr;
		}
	}
	if(_gridMode == KaiserBessel)
	{
		correctImageForKernel<false>(dataArray[0]);
	}
}

void LayeredImager::initializeSqrtLMLookupTable()
{
	_sqrtLMLookupTable.resize(_width * _height);
	std::vector<double>::iterator iter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - y) + _height / 2;
		if(ySrc >= _height) ySrc -= _height;
		double m = ((double) ySrc-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xSrc = x + _width / 2;
			if(xSrc >= _width) xSrc -= _width;
			double l = ((_width/2)-(double) xSrc) * _pixelSizeX + _phaseCentreDL;
			
			if(l*l + m*m < 1.0)
				*iter = sqrt(1.0 - l*l - m*m) - 1.0;
			else
				*iter = 0.0;
			++iter;
		}
	}
}

template<bool IsComplex>
void LayeredImager::projectOnImageAndCorrect(const std::complex<double> *source, double w, size_t threadIndex)
{
	double *dataReal = _imageData[threadIndex], *dataImaginary;
	if(IsComplex)
		dataImaginary = _imageDataImaginary[threadIndex];
	
	const double twoPiW = -2.0 * M_PI * w;
	std::vector<double>::const_iterator sqrtLMIter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - y) + _height / 2;
		if(ySrc >= _height) ySrc -= _height;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xSrc = x + _width / 2;
			if(xSrc >= _width) xSrc -= _width;
			
			double rad = twoPiW * *sqrtLMIter;
			double s, c;
			sincos(rad, &s, &c);
			/*std::complex<double> val = std::complex<double>(
				source->real() * c - source->imag() * s,
				source->real() * s + source->imag() * c
			);*/
			dataReal[xSrc + ySrc*_width] += source->real()*c - source->imag()*s;
			if(IsComplex)
			{
				if(_imageConjugatePart)
					dataImaginary[xSrc + ySrc*_width] += -source->real()*s + source->imag()*c;
				else
					dataImaginary[xSrc + ySrc*_width] += source->real()*s + source->imag()*c;
			}
			
			++source;
			++sqrtLMIter;
		}
	}
}

void LayeredImager::initializeSqrtLMLookupTableForSampling()
{
	_sqrtLMLookupTable.resize(_width * _height);
	std::vector<double>::iterator iter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		//size_t yDest = (_height - y) + _height / 2;
		size_t yDest = y + _height / 2;
		if(yDest >= _height) yDest -= _height;
		double m = ((double) yDest-(_height/2)) * _pixelSizeY + _phaseCentreDM;
		
		for(size_t x=0;x!=_width;++x)
		{
			//size_t xDest = x + _width / 2;
			size_t xDest = (_width - x) + _width / 2;
			if(xDest >= _width) xDest -= _width;
			double l = ((_width/2)-(double) xDest) * _pixelSizeX + _phaseCentreDL;
			
			if(l*l + m*m < 1.0)
				*iter = sqrt(1.0 - l*l - m*m) - 1.0;
			else
				*iter = 0.0;
			++iter;
		}
	}
}

template<bool IsComplex>
void LayeredImager::copyImageToLayerAndInverseCorrect(std::complex<double> *dest, double w)
{
	double *dataReal = _imageData[0], *dataImaginary;
	if(IsComplex)
		dataImaginary = _imageDataImaginary[0];
	
	const double twoPiW = 2.0 * M_PI * w;
	std::vector<double>::const_iterator sqrtLMIter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		// The fact that yDest is different than ySrc as above, is because of the way
		// fftw expects the data to be ordered.
		//size_t ySrc = (_height - y) + _height / 2;
		//if(ySrc >= _height) ySrc -= _height;
		size_t yDest = y + _height / 2;
		if(yDest >= _height) yDest -= _height;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xDest = (_width - x) + _width / 2;
			if(xDest >= _width) xDest -= _width;
			//size_t xSrc = x + _width / 2;
			//if(xSrc >= _width) xSrc -= _width;
			
			double rad = twoPiW * *sqrtLMIter;
			double s, c;
			sincos(rad, &s, &c);
			double realVal = dataReal[xDest + yDest*_width];
			if(IsComplex)
			{
				double imagVal = dataImaginary[xDest + yDest*_width];
				*dest = std::complex<double>(realVal*c + imagVal*s, imagVal*c - realVal*s);
			}
			else
				*dest = std::complex<double>(realVal*c, -realVal*s);
			
			++dest;
			++sqrtLMIter;
		}
	}
}

void LayeredImager::ReplaceRealImageBuffer(double* newBuffer)
{
	_imageBufferAllocator->Free(_imageData[0]);
	_imageData[0] = newBuffer;
}

void LayeredImager::ReplaceImaginaryImageBuffer(double* newBuffer)
{
	_imageBufferAllocator->Free(_imageDataImaginary[0]);
	_imageDataImaginary[0] = newBuffer;
}

