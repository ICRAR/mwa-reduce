#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <boost/thread/condition.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/function.hpp>

#include <unistd.h>
#include <queue>

class ThreadPool
{
public:
	ThreadPool() : _activeThreads(0), _finish(false)
	{
		init_threads(sysconf(_SC_NPROCESSORS_ONLN));
	}
	
	~ThreadPool()
	{
		boost::mutex::scoped_lock lock(_mutex);
		_finish = true;
		_change.notify_all();
		lock.unlock();
		for(size_t i=0; i!=_threads.size(); ++i)
		{
			_threads[i]->join();
		}
	}
	
	template<typename func>
	void queue(func f)
	{
		boost::mutex::scoped_lock lock(_mutex);
		_queuedTasks.push(f);
		_change.notify_all();
	}
	
	void wait_for_all_tasks()
	{
		boost::mutex::scoped_lock lock(_mutex);
		while(_activeThreads!=0 || !_queuedTasks.empty())
			_change.wait(lock);
	}
	
private:
	void init_threads(size_t n)
	{
		_threads.resize(n);
		for(size_t i=0; i!=n; ++i)
			_threads[i] = new boost::thread(
				std::mem_fun(&ThreadPool::thread_function), this);
	}
	
	void thread_function()
	{
		boost::mutex::scoped_lock lock(_mutex);
		while(!_finish)
		{
			while(_queuedTasks.empty() && !_finish)
				_change.wait(lock);
			
			if(!_finish)
			{
				++_activeThreads;
				boost::function<void()> f = _queuedTasks.front();
				lock.unlock();
				f();
				lock.lock();
				--_activeThreads;
			}
		}
	}
	
	std::vector<boost::thread*> _threads;
	std::size_t _activeThreads;
	std::queue<boost::function<void()> > _queuedTasks;
	bool _finish;
	
	boost::mutex _mutex;
	boost::condition _change;
};

#endif
