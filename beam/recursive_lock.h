#ifndef RECURSIVE_LOCK_H
#define RECURSIVE_LOCK_H

#include <cstring>
#include <mutex>
#include <system_error>

template<typename Mutex>
class recursive_lock
{
public:
	recursive_lock() :
		_mutex(nullptr),
		_nLocks(0)
	{ }
	
	recursive_lock(recursive_lock&& other) :
		_mutex(other._mutex),
		_nLocks(other._nLocks)
	{
		other._mutex = nullptr;
		other._nLocks = 0;
	}
	
	recursive_lock(Mutex& mutex) :
		_mutex(&mutex),
		_nLocks(1)
	{
		_mutex->lock();
	}
	
	recursive_lock(Mutex& mutex, std::defer_lock_t) noexcept :
		_mutex(&mutex),
		_nLocks(0)
	{ }
	
	recursive_lock(Mutex& mutex, std::try_to_lock_t) :
		_mutex(&mutex),
		_nLocks(0)
	{
		try_lock();
	}
	
	recursive_lock(Mutex& mutex, std::adopt_lock_t) noexcept :
		_mutex(&mutex),
		_nLocks(1)
	{	}
	
	~recursive_lock()
	{
		if(_nLocks != 0)
			_mutex->unlock();
	}
	
	recursive_lock& operator=(const recursive_lock& other)
	{
		if(_nLocks != 0)
			_mutex->unlock();
		_mutex = other._mutex;
		_nLocks = other.nLocks;
		other._mutex = nullptr;
		other._nLocks = 0;
	}
	
	void lock()
	{
		if(_mutex == nullptr)
			throw std::system_error(std::make_error_code(std::errc::operation_not_permitted), "lock() called on recursive_lock without mutex");
		if(_nLocks == 0)
			_mutex->lock();
		++_nLocks;
	}
	
	bool try_lock()
	{
		if(_mutex == nullptr)
			throw std::system_error(std::make_error_code(std::errc::operation_not_permitted), "lock() called on recursive_lock without mutex");
		if(_nLocks == 0)
		{
			if(_mutex->try_lock())
			{
				++_nLocks;
				return true;
			}
			else {
				return false;
			}
		}
		else {
			++_nLocks;
			return true;
		}
	}
	
	void unlock()
	{
		_nLocks--;
		if(_mutex == nullptr)
			throw std::system_error(std::make_error_code(std::errc::operation_not_permitted), "unlock() called on recursive_lock without mutex");
		else if(_nLocks == 0)
			_mutex->unlock();
	}
	
	void swap(recursive_lock& other) noexcept
	{
		std::swap(_mutex, other._mutex);
		std::swap(_nLocks, other._nLocks);
	}
	
	bool owns_lock() const noexcept
	{
		return _nLocks != 0;
	}
	
	Mutex* release() noexcept {
		Mutex* m = _mutex;
		_mutex = nullptr;
		_nLocks = 0;
		return m;
	}
	
	Mutex* mutex() const noexcept { return _mutex; }
	
	operator bool() const noexcept { return owns_lock(); }
	
private:
	Mutex* _mutex;
	size_t _nLocks;
};

#endif
