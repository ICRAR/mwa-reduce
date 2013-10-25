#include "lane_11.h"
#include "stopwatch.h"

#include <iostream>
#include <thread>

using namespace std;

const size_t NITER = 5000000;
const size_t N_AT_A_TIME = 1000;

void producer1(ao::lane<int>* lane)
{
	for(int a=0; a!=int(NITER); ++a)
		lane->write(a);
	lane->write_end();
}

void consumer1(ao::lane<int>* lane)
{
	int a, b = 0;
	while(lane->read(a))
	{
		if(a != b)
			std::cout << a << " != " << b << '\n';
		++b;
	}
	if(b != int(NITER))
		std::cout << "Incorrect nr of samples: got " << b << '\n';
}

void producerN(ao::lane<int>* lane, size_t n)
{
	std::vector<int> vals(n);
	int i = 0;
	for(int a=0; a!=int(NITER/n); ++a)
	{
		for(int b=0; b!=int(n); ++b)
		{
			vals[b] = i;
		  ++i;
		}
		lane->write(vals.data(), n);
	}
	for(size_t b=0; b!=NITER%n; ++b)
	{
		vals[b] = i;
		++i;
	}
	lane->write(vals.data(), NITER%n);
	lane->write_end();
}

void consumerN(ao::lane<int>* lane, size_t n)
{
	std::vector<int> vals(n);
	int b = 0;
	size_t nread = lane->read(vals.data(), n);
	while(nread != 0)
	{
		for(size_t i=0; i!=nread; ++i)
		{
			if(vals[i] != b)
				std::cout << vals[i] << " != " << b << '\n';
			++b;
		}
		nread = lane->read(vals.data(), n);
	}
	if(b != int(NITER))
		std::cout << "Incorrect nr of samples: got " << b << '\n';
}

int main(int argc, char* argv[])
{
	std::cout << "One at a time read & write: " << std::flush;
	ao::lane<int> l(N_AT_A_TIME*2);
	Stopwatch watch(true);
	thread producerThread1(&producer1, &l);
	thread consumerThread1(&consumer1, &l);
	producerThread1.join();
	consumerThread1.join();
	std::cout << watch.ToString() << " (" << (NITER/watch.Seconds()) << " read-writes/sec)\n";
	
	std::cout << "N at a time write, one at a time read: " << std::flush;
	l.clear();
	watch.Reset();
	watch.Start();
	thread producerThread2(&producerN, &l, N_AT_A_TIME);
	thread consumerThread2(&consumer1, &l);
	producerThread2.join();
	consumerThread2.join();
	std::cout << watch.ToString() << " (" << (NITER/watch.Seconds()) << " read-writes/sec)\n";
	
	std::cout << "N at a time read, one at a time write: " << std::flush;
	l.clear();
	watch.Reset();
	watch.Start();
	thread producerThread3(&producer1, &l);
	thread consumerThread3(&consumerN, &l, N_AT_A_TIME);
	producerThread3.join();
	consumerThread3.join();
	std::cout << watch.ToString() << " (" << (NITER/watch.Seconds()) << " read-writes/sec)\n";
	
	std::cout << "N at a time read & write: " << std::flush;
	l.clear();
	watch.Reset();
	watch.Start();
	thread producerThread4(&producerN, &l, N_AT_A_TIME);
	thread consumerThread4(&consumerN, &l, N_AT_A_TIME);
	producerThread4.join();
	consumerThread4.join();
	std::cout << watch.ToString() << " (" << (NITER/watch.Seconds()) << " read-writes/sec)\n";
	
	std::cout << "N at a time read & write (write from this thread): " << std::flush;
	l.clear();
	watch.Reset();
	watch.Start();
	thread consumerThread5(&consumerN, &l, N_AT_A_TIME);
	producerN(&l, N_AT_A_TIME);
	consumerThread5.join();
	std::cout << watch.ToString() << " (" << (NITER/watch.Seconds()) << " read-writes/sec)\n";
	
	size_t n = 1;
	while(n < 100000)
	{
		l.clear();
		watch.Reset();
		watch.Start();
		thread consumerThread6(&consumerN, &l, n);
		producerN(&l, n);
		consumerThread6.join();
		std::cout << n << '\t' << (NITER/watch.Seconds()) << '\n';
		
		n *= 2;
	}
	
	return 0;
}
