#include "layeredimager.h"

#include <iostream>

#include <fftw3.h>

LayeredImager::LayeredImager(size_t width, size_t height, double pixelScale) :
	_width(width),
	_height(height),
	_pixelScale(pixelScale)
{
	size_t imgSize = _height * _width;
	posix_memalign(reinterpret_cast<void**>(&_imageData), sizeof(double)*2, imgSize * sizeof(double));
	memset(_imageData, 0, imgSize * sizeof(double));
	
	initializeSqrtLMLookupTable();
}

LayeredImager::~LayeredImager()
{
	free(_imageData);
}

void LayeredImager::PrepareForObservation(size_t nWLayers, size_t maxMem, double minW, double maxW, const BandData &bandData)
{
	_minW = minW;
	_maxW = maxW;
	_nWLayers = nWLayers;
	
	size_t maxNWLayersPerPass = size_t((double) maxMem / (_width * _height * sizeof(double) * 2));
	_nPasses = (nWLayers+maxNWLayersPerPass-1)/maxNWLayersPerPass;
	if(_nPasses == 0) _nPasses = 1;
	std::cout << "Will process " << (_nWLayers / _nPasses) << " w-layers per pass.\n";
	
	_curLayerRangeIndex = 0;
	_bandData = &bandData;
}

void LayeredImager::StartPass(size_t passIndex)
{
	_curLayerRangeIndex = passIndex;
	size_t nLayersInPass = layerRangeStart(passIndex+1) - layerRangeStart(passIndex);
	_layeredUVData.resize(nLayersInPass);
	for(size_t l=0; l!=nLayersInPass; ++l)
	{
		_layeredUVData[l].assign(_width * _height, std::complex<double>());
	}
}

void LayeredImager::FinishPass()
{
	size_t imgSize = _width * _height;
	std::complex<double> *fftwIn = reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	std::complex<double> *fftwOut = reinterpret_cast<std::complex<double>*>(fftw_malloc(imgSize * sizeof(double) * 2));
	
	fftw_plan plan =
		fftw_plan_dft_2d(_width, _height,
			reinterpret_cast<fftw_complex*>(fftwIn), reinterpret_cast<fftw_complex*>(fftwOut),
			FFTW_BACKWARD, FFTW_ESTIMATE);

	size_t layerOffset = layerRangeStart(_curLayerRangeIndex);
	size_t nPlanes = layerRangeStart(_curLayerRangeIndex+1) - layerOffset;
	for(size_t plane=0; plane!=nPlanes; ++plane)
	{
		// Fourier transform the layer
		std::vector<std::complex<double>> &uvData = _layeredUVData[plane];
		memcpy(fftwIn, &uvData[0], imgSize * sizeof(double) * 2);
		fftw_execute(plan);
		
		// Add layer to full image
		projectOnImageAndCorrect(fftwOut, LayerToW(plane + layerOffset));
	}
	
	fftw_free(fftwIn);
	fftw_free(fftwOut);
}

void LayeredImager::AddData(const std::complex<float>* data, double uInM, double vInM, double wInM)
{
 	const size_t
		layerOffset = layerRangeStart(_curLayerRangeIndex),
		layerRangeEnd = layerRangeStart(_curLayerRangeIndex+1);
	for(size_t ch=0; ch!=_bandData->ChannelCount(); ++ch)
	{
		double
			wavelength = _bandData->ChannelWavelength(ch),
			u = uInM / wavelength,
			v = vInM / wavelength,
			w = wInM / wavelength;
		std::complex<float> sample = data[ch];
		if(w < 0.0)
		{
			u = -u;
			v = -v;
			w = -w;
			sample = std::conj(sample);
		}
		size_t
			wLayer = WToLayer(w);
		if(wLayer >= layerOffset && wLayer < layerRangeEnd)
		{
			size_t layerIndex = wLayer - layerOffset;
			std::vector<std::complex<double>> &uvData = _layeredUVData[layerIndex];
			int
				x = int(round(u * _pixelScale * _width)),
				y = int(round(v * _pixelScale * _height));
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

void LayeredImager::FinalizeImage()
{
	_layeredUVData.clear();
	std::vector<double> image(_width * _height);
	for(size_t y=0;y!=_height;++y)
	{
		for(size_t x=0;x!=_width;++x)
		{
			image[x + y*_width] *= 1.0; //TODO
		}
	}
}

void LayeredImager::initializeSqrtLMLookupTable()
{
	_sqrtLMLookupTable.resize(_width * _height);
	std::vector<double>::iterator iter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - 1 - y) + _height / 2;
		if(ySrc >= _height) ySrc -= _height;
		double m = ((double) ySrc-(_height/2)) * _pixelScale;
		
		for(size_t x=0;x!=_width;++x)
		{
			size_t xSrc = x + _width / 2;
			if(xSrc >= _width) xSrc -= _width;
			
			double l = ((double) xSrc-(_width/2)) * _pixelScale;
			*iter = (sqrt(1.0 - l*l - m*m)-1.0);
			iter++;
		}
	}
}

void LayeredImager::projectOnImageAndCorrect(const std::complex<double> *source, double w)
{
	const double twoPiW = -2.0 * M_PI * w;
	std::vector<double>::const_iterator sqrtLMIter = _sqrtLMLookupTable.begin();
	for(size_t y=0;y!=_height;++y)
	{
		size_t ySrc = (_height - 1 - y) + _height / 2;
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
			_imageData[xSrc + ySrc*_width] += source->real()*c - source->imag()*s;
			
			++source;
			++sqrtLMIter;
		}
	}
}
