#ifndef MSIOSYSTEM_H
#define MSIOSYSTEM_H

#include <casacore/casa/OS/HostInfo.h>

#include <stdio.h>
#include <unistd.h>
#include <sched.h>
#include <string.h>

#include <stdexcept>
#include <string>

class System
{
	public:
		static long TotalMemory()
		{
			return casacore::HostInfo::memoryTotal()*1024;
		}
		
		static unsigned ProcessorCount()
		{
#ifdef __APPLE__
			return sysconf(_SC_NPROCESSORS_ONLN);
#else
			cpu_set_t cs;
			CPU_ZERO(&cs);
			sched_getaffinity(0, sizeof cs , &cs);

			int count = 0;
			for (int i = 0; i < CPU_SETSIZE; i++)
			{
				if (CPU_ISSET(i, &cs))
				count++;
			}

			return count;
#endif
		}
		
	private:
		static char* handle_strreturn(int value)
		{
			if(value != 0)
				throw std::runtime_error("strerror_r() reported an error");
			return 0;
		}
		static char* handle_strreturn(char* value)
		{
			return value;
		}
	public:
		static std::string StrError(int errnum)
		{
			// Because strerror_r() has different return values on different platforms,
			// two overloads of handle_strerror are used to make this compile and work
			// in either case of int or char*.
			char buffer[1024];
			char* ret = handle_strreturn(strerror_r(errnum, buffer, 1024));
			if(ret == 0)
				return std::string(buffer);
			else
				return std::string(ret);
		}
};

#endif //MSIOSYSTEM_H
