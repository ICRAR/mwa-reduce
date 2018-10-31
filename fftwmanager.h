#ifndef FFTW_MULTI_THREAD_ENABLER_H
#define FFTW_MULTI_THREAD_ENABLER_H

#include <cstring>

#include <mutex>

/**
 * Used to initialize and enable fftw's multithreading. While
 * an instance of this class exists, multi-threading is enabled. When the
 * class is destructed, multi-threading is again disabled.
 * 
 * To make the FFTs in for example the @ref FFTConvolver or @ref FFTResampler
 * multi-threaded, it is enough to construct an instance of this class and keep
 * it until done.
 * For example:
 * 
 * @code
 * void doMultiThreadedFFTActions
 * {
 *   FFTWManager fftwMT;
 *   FFTWManager::ThreadingScope fftwScope(fftwMT);
 * 
 *   FFTResampler resampler(..);
 *   ...
 * }
 * @endcode
 */
class FFTWManager
{
public:
	class ThreadingScope
	{
	public:
		ThreadingScope(FFTWManager& manager) :
			_manager(manager)
		{
			manager.IncreaseMultiThreadUsers();
		}
		~ThreadingScope()
		{
			_manager.DecreaseMultiThreadUsers();
		}
	private:
		FFTWManager& _manager;
		ThreadingScope(const ThreadingScope&) = delete;
		ThreadingScope* operator=(const ThreadingScope&) = delete;
	};
	
	/**
	 * Constructor that sets FFTW to use multiple threads.
	 * This will set FFTW to use as many threads as there are cores in the 
	 * system.
	 */
	explicit FFTWManager(bool verbose = false);
	
	/**
	 * Constructor that sets FFTW to use multiple threads.
	 * This will set FFTW to use the given number of threads.
	 */
	explicit FFTWManager(size_t nThreads, bool verbose = false);
	
	~FFTWManager();
	
	void IncreaseMultiThreadUsers()
	{
		std::lock_guard<std::mutex> lock(_mutex);
		if(_multiThreadEnabledDepth == 0)
			activateMultipleThreads();
		++_multiThreadEnabledDepth;
	}
	
	void DecreaseMultiThreadUsers()
	{
		std::lock_guard<std::mutex> lock(_mutex);
		--_multiThreadEnabledDepth;
		if(_multiThreadEnabledDepth == 0)
			endMultipleThreads();
	}
	
	std::mutex& Mutex() { return _mutex; }
	
private:
	void activateMultipleThreads();
	void endMultipleThreads();
	
	int _multiThreadEnabledDepth;
	bool _verbose;
	size_t _nThreads;
	std::mutex _mutex;
};

#endif
