#ifndef RECURSIVE_LOCK_H
#define RECURSIVE_LOCK_H

#include <cstring>

template<typename Mutex>
class recursive_lock
{
public:
	recursive_lock(Mutex& mutex, bool acquire_lock = true) :
		_mutex(&mutex)
	{
		if(acquire_lock)
		{
			_mutex->lock();
			_nLocks = 1;
		}
		else {
			_nLocks = 0;
		}
	}
	
	~recursive_lock()
	{
		if(_nLocks != 0)
			_mutex->unlock();
	}
	
	void lock()
	{
		if(_nLocks == 0)
			_mutex->lock();
		++_nLocks;
	}
	
	void unlock()
	{
		_nLocks--;
		if(_nLocks == 0)
			_mutex->unlock();
	}
	
private:
	Mutex* _mutex;
	size_t _nLocks;
};

#endif
