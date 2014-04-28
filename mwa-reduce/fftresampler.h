#ifndef FFT_RESAMPLE_H
#define FFT_RESAMPLE_H

#include "lane.h"

#include <thread>
#include <vector>

#include <fftw3.h>

class FFTResampler
{
private:
	struct Task
	{
		double *input, *output;
	};
	
public:
	FFTResampler(size_t inWidth, size_t inHeight, size_t outWidth, size_t outHeight, size_t cpuCount, bool verbose=true);
	
	~FFTResampler()
	{
		Finish();
	}
	
	void AddTask(double* input, double* output)
	{
		Task task;
		task.input = input;
		task.output = output;
		_tasks.write(task);
	}
	
	void Start()
	{
		for(size_t i=0; i!=_tasks.capacity(); ++i)
		{
			_threads.push_back(std::thread(&FFTResampler::runThread, this));
		}
	}
	
	void Finish()
	{
		_tasks.write_end();
		for(std::vector<std::thread>::iterator i=_threads.begin(); i!=_threads.end(); ++i)
			i->join();
		_threads.clear();
		_tasks.clear();
	}
	
	void RunSingle(double* input, double* output)
	{
		AddTask(input, output);
		_tasks.write_end();
		runThread();
		_tasks.clear();
	}
	
private:
	void runThread();
	
	size_t _inputWidth, _inputHeight;
	size_t _outputWidth, _outputHeight;
	size_t _fftWidth, _fftHeight;
	
	fftw_plan _inToFPlan, _fToOutPlan;
	
	ao::lane<Task> _tasks;
	std::vector<std::thread> _threads;
	bool _verbose;
};

#endif
