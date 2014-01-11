#ifndef SOLUTION_FILE_H
#define SOLUTION_FILE_H

#include <cstring>
#include <string>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <mutex>

#include <stdint.h>

class IonSolutionFile
{
 public:
	 enum IonSolutionType { GainSolution, DlSolution, DmSolution };
	 
	 struct Solution
	 {
		 double gain, dl, dm;
	 };
	 
  IonSolutionFile() : _outputStream(0), _inputStream(0)
  {
    strcpy(_header.intro, "MWAOCAL");
    _header.fileType = 1; // Multi-directional gain,posl,posm solutions
    _header.structureType = 0; // ordered gain/dl/dm, polarization, channel, antenna, time
    _header.intervalCount = 1;
		_header.antennaCount = 1;
		_header.channelCount = 0;
		_header.polarizationCount = 1;
		_header.directionCount = 0;
		_header.parameterCount = 3;
		_header.startTime = 0.0;
		_header.endTime = 0.0;
  }

  ~IonSolutionFile() {
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

  size_t DirectionCount() const { return _header.directionCount; }
  void SetDirectionCount(size_t directionCount) {
		_header.directionCount = directionCount;
	}

  void OpenForWriting(const char *filename)
  {
		delete _outputStream;
		_outputStream = new std::ofstream(filename);    
		_outputStream->write(reinterpret_cast<const char*>(&_header), sizeof(_header));
  }
  
	void OpenForReading(const char *filename)
	{
		delete _inputStream;
		_inputStream = new std::ifstream(filename);
		if(_inputStream->bad())
			throw std::runtime_error("Error reading input ionospheric solutions file");
		_inputStream->read(reinterpret_cast<char*>(&_header), sizeof(_header));
	}

	double ReadSolution(IonSolutionType type, size_t interval, size_t channel, size_t polarization, size_t direction)
	{
		Solution solution;
		ReadSolution(solution, interval, channel, polarization, direction);
		switch(type)
		{
			default:
			case GainSolution: return solution.gain;
			case DlSolution: return solution.dl;
			case DmSolution: return solution.dm;
		}
	}
	
  void ReadSolution(Solution& solution, size_t interval, size_t channel, size_t polarization, size_t direction) {
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelCount + channel) * _header.polarizationCount + polarization) * _header.directionCount + direction;
		size_t offset = sizeof(_header);
		_inputStream->seekg(offset + sizeof(Solution) * index, std::ios::beg);
		_inputStream->read(reinterpret_cast<char*>(&solution), sizeof(Solution));
  }

  void WriteSolution(const Solution& solution, size_t interval, size_t channel, size_t polarization, size_t direction)
  {
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelCount + channel) * _header.polarizationCount + polarization) * _header.directionCount + direction;
		size_t offset = sizeof(_header);
		_outputStream->seekp(offset + sizeof(Solution) * index, std::ios::beg);
		_outputStream->write(reinterpret_cast<const char*>(&solution), sizeof(Solution));
  }

  void WriteChannelBlock(const Solution* solutions, size_t interval, size_t channel, size_t polarization)
  {
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelCount + channel) * _header.polarizationCount + polarization) * _header.directionCount;
		size_t offset = sizeof(_header);
		_outputStream->seekp(offset + sizeof(Solution) * index, std::ios::beg);
		_outputStream->write(reinterpret_cast<const char*>(solutions), sizeof(Solution) * _header.directionCount);
  }

 private:
  struct {
    char intro[8];
    uint32_t fileType;
    uint32_t structureType;
    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;
		uint32_t directionCount, parameterCount;
		double startTime, endTime;
  } _header;
  std::ofstream *_outputStream;
  std::ifstream *_inputStream;
	std::mutex _mutex;
};

#endif
