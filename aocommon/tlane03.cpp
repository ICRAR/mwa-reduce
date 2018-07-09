#include "lane_03.h"
#include "stopwatch.h"

#include <iostream>
#include <thread>

using namespace std;

const size_t NITER = 5000000;

void producer(ao::lane<int>* lane)
{
	for(int a=0; a!=int(NITER); ++a)
		lane->write(a);
	lane->write_end();
}

void consumer(ao::lane<int>* lane)
{
	int a;
	while(lane->read(a))
	{ }
}

int main(int argc, char* argv[])
{
	ao::lane<int> l(1000);
	
	Stopwatch watch(true);
	thread producerThread(&producer, &l);
	thread consumerThread(&consumer, &l);
	producerThread.join();
	consumerThread.join();
	std::cout << "Time=" << watch.ToString() << '\n';
	std::cout << (NITER/watch.Seconds()) << " read-writes/sec\n";
	return 0;
}
