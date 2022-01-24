#ifndef SOLUTION_FLAG_FILE_H
#define SOLUTION_FLAG_FILE_H

#include <aocommon/uvector.h>

#include <string>
#include <fstream>

class SolutionFlagFile {
 public:
  SolutionFlagFile(const std::string& filename)
      : _intervalCount(0), _antennaCount(0), _channelCount(0) {
    std::ifstream str(filename);
    if (!str.good())
      throw std::runtime_error("Could not open specified solution flag file");
    std::string line;
    std::getline(str, line);
    _intervalCount = 1;
    while (str.good()) {
      std::string line;
      std::getline(str, line);
      size_t index = 0;
      if (!line.empty()) {
        size_t cCount = (line.size() + 1) / 2;
        if (_channelCount == 0)
          _channelCount = cCount;
        else if (_channelCount != cCount)
          throw std::runtime_error(
              "Format of specified solution flag file is not consistent");
        _flags.push_back(_channelCount, false);
        for (size_t ch = 0; ch != _channelCount; ++ch) {
          bool f;
          if (line[ch * 2] == '0')
            f = false;
          else if (line[ch * 2] == '1')
            f = true;
          else
            throw std::runtime_error("Expecting 0 or 1 in flag file");
          if (ch != 0 && line[ch * 2 - 1] != ',')
            throw std::runtime_error(
                "Expecting commas between 0/1s in flag file");
          _flags[index + ch] = f;
        }
        index += _channelCount;
        ++_antennaCount;
      }
    }
  }

  size_t IntervalCount() const { return _intervalCount; }
  size_t AntennaCount() const { return _antennaCount; }
  size_t ChannelCount() const { return _channelCount; }

  bool IsFlagged(size_t intervalIndex, size_t antennaIndex,
                 size_t channelIndex) const {
    return _flags[antennaIndex * _channelCount + channelIndex];
  }

 private:
  aocommon::UVector<bool> _flags;
  size_t _intervalCount, _antennaCount, _channelCount;
};

#endif
