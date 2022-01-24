#include "progressbar.h"

#include <iostream>

ProgressBar::ProgressBar(const std::string& taskDescription)
    : _taskDescription(taskDescription), _displayedDots(-1) {}

ProgressBar::~ProgressBar() { SetProgress(1, 1); }

ProgressBar& ProgressBar::operator=(ProgressBar&& rhs) {
  SetProgress(1, 1);
  _displayedDots = rhs._displayedDots;
  _taskDescription = std::move(rhs._taskDescription);
  rhs._displayedDots = 50;
  return *this;
}

void ProgressBar::SetProgress(size_t taskIndex, size_t taskCount) {
  if (_displayedDots == -1) {
    std::cout << _taskDescription << ":";
    if (_taskDescription.size() < 40)
      std::cout << " 0%";
    else
      std::cout << "\n 0%";
    _displayedDots = 0;
    std::cout.flush();
  }
  int progress = (taskIndex * 100 / taskCount);
  int dots = progress / 2;

  if (dots > _displayedDots) {
    while (dots != _displayedDots) {
      ++_displayedDots;
      if (_displayedDots % 5 == 0)
        std::cout << ((_displayedDots / 5) * 10) << '%';
      else
        std::cout << '.';
    }
    if (progress == 100) std::cout << '\n';
    std::cout.flush();
  }
}
