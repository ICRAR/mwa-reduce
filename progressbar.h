#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <string>

class ProgressBar
{
	public:
		explicit ProgressBar(const std::string &taskDescription);
		~ProgressBar();
	
		void SetProgress(size_t taskIndex, size_t taskCount);
		
		ProgressBar& operator=(ProgressBar&& rhs);
		
	private:
		std::string _taskDescription;
		int _displayedDots;
};

#endif
