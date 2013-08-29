#ifndef SOLUTION_FILE_H
#define SOLUTION_FILE_H

#include <cstring>
#include <string>
#include <fstream>
#include <complex>
#include <stdexcept>

#include <stdint.h>

class SolutionFile
{
 public:
  SolutionFile() : _outputStream(0), _inputStream(0)
  {
    strcpy(_header.intro, "MWAOCAL");
    _header.fileType = 0; // Complex jones solutions
    _header.structureType = 0; // ordered real/imag, polarization, channel, antenna, time
    _header.intervalCount = 1;
  }

  ~SolutionFile() {
    delete _outputStream;
    delete _inputStream;
  }

  size_t AntennaCount() const { return _header.antennaCount; }
  void SetAntennaCount(size_t antennaCount) {
    _header.antennaCount = antennaCount;
  }

  size_t ChannelCount() const { return _header.channelCount; }
  void SetChannelCount(size_t channelCount) {
    _header.channelCount = channelCount;
  }

  size_t PolarizationCount() const { return _header.polarizationCount; }
  void SetPolarizationCount(size_t polarizationCount) {
    _header.polarizationCount = polarizationCount;
  }
  
  size_t IntervalCount() const { return _header.intervalCount; }
  void SetIntervalCount(size_t intervalCount) {
		_header.intervalCount = intervalCount;
	}

  void OpenForWriting(const char *filename)
  {
    delete _outputStream;
    _outputStream = new std::ofstream(filename);    
    _outputStream->write(reinterpret_cast<const char*>(&_header), sizeof(_header));
    double timeStart = 0.0, timeEnd = 0.0;
    _outputStream->write(reinterpret_cast<const char*>(&timeStart), sizeof(timeStart));
    _outputStream->write(reinterpret_cast<const char*>(&timeEnd), sizeof(timeEnd)); 
  }

  void OpenForReading(const char *filename)
  {
    delete _inputStream;
    _inputStream = new std::ifstream(filename);
    if(_inputStream->bad())
      throw std::runtime_error("Error reading input solutions file");
    _inputStream->read(reinterpret_cast<char*>(&_header), sizeof(_header));
    double timeStart, timeEnd;
    _inputStream->read(reinterpret_cast<char*>(&timeStart), sizeof(timeStart));
    _inputStream->read(reinterpret_cast<char*>(&timeEnd), sizeof(timeEnd)); 
  }

  std::complex<double> ReadNextSolution() {
    std::complex<double> val;
    _inputStream->read(reinterpret_cast<char*>(&val), sizeof(val));
    return val;
  }

  void WriteNextSolution(const std::complex<double> &val)
  {
    _outputStream->write(reinterpret_cast<const char*>(&val), sizeof(val));
  }

  void WriteSolution(const std::complex<double> &val, size_t interval, size_t antenna, size_t channel, size_t polarization)
  {
		size_t index = ((interval * _header.intervalCount + antenna) * _header.channelCount + channel) * _header.polarizationCount + polarization;
		size_t offset = sizeof(_header) + sizeof(double)*2;
		_outputStream->seekp(offset + sizeof(val) * index, std::ios::beg);
    _outputStream->write(reinterpret_cast<const char*>(&val), sizeof(val));
  }

 private:
  struct {
    char intro[8];
    uint32_t fileType;
    uint32_t structureType;
    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;
  } _header;
  std::ofstream *_outputStream;
  std::ifstream *_inputStream;
};

#endif
