#ifndef SYNCHRONIZED_MS_H
#define SYNCHRONIZED_MS_H

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <condition_variable>
#include <mutex>
#include <set>
#include <string>

class SynchronizedMS
{
public:
	SynchronizedMS()
	{ }
	
	SynchronizedMS(const std::string& filename) :
		_lock(std::make_shared<MSLock>(filename))
	{ }
	
	void Reset()
	{
		_lock.reset();
	}
	
	casacore::MeasurementSet* operator->() { return &_lock->MS(); }
	
	casacore::MeasurementSet& operator*() { return _lock->MS(); }
	
	casacore::MeasurementSet& MS() { return _lock->MS(); }
	
private:
	
	class MSLock {
	public:
		MSLock(const std::string& filename) : _filename(filename)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			while(_openFiles.count(filename))
				_condition.wait(lock);
			_openFiles.insert(filename);
			lock.unlock();
			
			_ms.reset(new casacore::MeasurementSet(_filename));
		}
		
		MSLock(const MSLock&) = delete;
		
		MSLock(MSLock&& other)
		{
			// the other must already have access, so the ms can be directly claimed.
			// (unless the other was an empty object, in which case the
			//  copy will also be an empty object)
			_filename = std::move(other._filename);
			_ms = std::move(other._ms);
			other._filename = std::string();
			other._ms.reset();
		}
		
		~MSLock()
		{
			release();
		}
		
		MSLock& operator=(const MSLock&) = delete;
		
		MSLock& operator=(MSLock&& rhs)
		{
			release();
			
			_filename = std::move(rhs._filename);
			_ms = std::move(rhs._ms);
			return *this;
		}
		
		casacore::MeasurementSet& MS() { return *_ms; }
		
	private:
		void release()
		{
			if(!_filename.empty())
			{
				_ms.reset();
				std::lock_guard<std::mutex> lock(_mutex);
				_openFiles.erase(_filename);
				_condition.notify_all();
			}
		}
		
		static std::set<std::string> _openFiles;
		static std::condition_variable _condition;
		static std::mutex _mutex;
		
		std::string _filename;
		std::unique_ptr<casacore::MeasurementSet> _ms;
	};
	
	std::shared_ptr<MSLock> _lock;
};

#endif
