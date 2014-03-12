#include "fftresampler.h"

#include <complex>
#include <iostream>

FFTResampler::FFTResampler(size_t inWidth, size_t inHeight, size_t outWidth, size_t outHeight, size_t cpuCount) :
	_inputWidth(inWidth), _inputHeight(inHeight),
	_outputWidth(outWidth), _outputHeight(outHeight),
	_fftWidth(std::max(inWidth, outWidth)), _fftHeight(std::max(inHeight, outHeight)),
	_tasks(cpuCount)
{
	double* inputData = reinterpret_cast<double*>(fftw_malloc(_fftWidth*_fftHeight * sizeof(double)));
	fftw_complex* fftData = reinterpret_cast<fftw_complex*>(fftw_malloc(_fftWidth*_fftHeight * sizeof(fftw_complex)));
	_inToFPlan =
		fftw_plan_dft_r2c_2d(_inputHeight, _inputWidth,
			inputData, fftData, FFTW_ESTIMATE);
	_fToOutPlan =
		fftw_plan_dft_c2r_2d(_outputHeight, _outputWidth,
			fftData, inputData, FFTW_ESTIMATE);
	fftw_free(fftData);
	fftw_free(inputData);
}

void FFTResampler::runThread()
{
	Task task;
	while(_tasks.read(task))
	{
		double *endPtr = task.input + _inputWidth*_inputHeight;
		for(double* i=task.input; i!=endPtr; ++i)
		{
			if(!std::isfinite(*i))
				*i = 0.0;
		}
		
		size_t fftInWidth = _inputWidth/2+1;
		std::complex<double>
			*fftData = reinterpret_cast<std::complex<double>*>(fftw_malloc(fftInWidth*_inputHeight*sizeof(std::complex<double>)));
		std::cout << "FFT " << _inputWidth << " x " << _inputHeight << " real -> complex...\n";
		fftw_execute_dft_r2c(_inToFPlan, task.input, reinterpret_cast<fftw_complex*>(fftData));
		
		size_t fftOutWidth = _outputWidth/2+1;
		// TODO this can be done without allocating more mem!
		std::complex<double>
			*newfftData = reinterpret_cast<std::complex<double>*>(fftw_malloc(fftOutWidth*_outputHeight*sizeof(std::complex<double>)));
		memset(newfftData, 0, fftOutWidth*_outputHeight*sizeof(std::complex<double>));
			
		size_t oldMidX = _inputWidth/2;
		size_t newMidX = _outputWidth/2;
		
		size_t minWidth = std::min(_inputWidth, _outputWidth);
		size_t minHeight = std::min(_inputHeight, _outputHeight);
		
		size_t minMidX = minWidth/2;
		size_t minMidY = minHeight/2;
		
		double factor = 1.0 / (minWidth*minHeight);
		
		for(size_t y=0; y!=minHeight; ++y)
		{
			size_t oldY = y-minMidY + _inputHeight;
			size_t newY = y-minMidY + _outputHeight;
			if(oldY >= _inputHeight) oldY -= _inputHeight;
			if(newY >= _outputHeight) newY -= _outputHeight;
			
			// The last dimension is stored half
			for(size_t x=0; x!=minMidX; ++x)
			{
				size_t oldX = x;
				size_t newX = x;
				size_t oldIndex = oldX + oldY * (oldMidX+1);
				size_t newIndex = newX + newY * (newMidX+1);
				
				newfftData[newIndex] = fftData[oldIndex] * factor;
				
				//if((x == 0 && newY == 0) || (x==minMidX-1 && y==minHeight-1))
				//	std::cout << newfftData[newIndex] << " (" << oldX << " , " << oldY << ") - (" << newX << " , " << newY << ")\n";
			}
			if(_inputWidth > _outputWidth)
			{
				size_t oldIndex = _inputWidth/2 + oldY * (oldMidX+1);
				size_t newIndex = _outputWidth/2 + newY * (newMidX+1);
				newfftData[newIndex] = fftData[oldIndex] * factor;
			}
		}
		
		fftw_free(fftData);
		
		std::cout << "FFT " << _outputWidth << " x " << _outputHeight << " complex -> real...\n";
		fftw_execute_dft_c2r(_fToOutPlan, reinterpret_cast<fftw_complex*>(newfftData), task.output);
		
		fftw_free(newfftData);
	}
}
