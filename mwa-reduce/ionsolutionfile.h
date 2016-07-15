#ifndef SOLUTION_FILE_H
#define SOLUTION_FILE_H

#include <cstring>
#include <string>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <mutex>
#include <cstdint>

#include "uvector.h"

#include "model/model.h"

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
		// 1: Multi-directional gain,posl,posm solutions, old format without cluster descs
		// 2: Same, but new format including cluster descs
    _header.fileType = 2;
    _header.structureType = 0; // ordered gain/dl/dm, polarization, channel, antenna, time
    _header.intervalCount = 1;
		_header.antennaCount = 1;
		_header.channelBlockCount = 0;
		_header.polarizationCount = 1;
		_header.directionCount = 0;
		_header.parameterCount = 3;
		_header.startTime = 0.0;
		_header.endTime = 0.0;
		_header.startFrequency = 0.0;
		_header.endFrequency = 0.0;
  }

  ~IonSolutionFile() {
    delete _outputStream;
    delete _inputStream;
  }

  size_t AntennaCount() const { return _header.antennaCount; }
  void SetAntennaCount(size_t antennaCount) {
    _header.antennaCount = antennaCount;
  }

  size_t ChannelBlockCount() const { return _header.channelBlockCount; }
  void SetChannelBlockCount(size_t channelBlockCount) {
    _header.channelBlockCount = channelBlockCount;
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
	
	double StartFrequency() const { return _header.startFrequency; }
	void SetStartFrequency(double startFrequency) { _header.startFrequency = startFrequency; }
	
	double EndFrequency() const { return _header.endFrequency; }
	void SetEndFrequency(double endFrequency) { _header.endFrequency = endFrequency; }

  void OpenForWriting(const char *filename)
  {
		delete _outputStream;
		_outputStream = new std::ofstream(filename);    
		_outputStream->write(reinterpret_cast<const char*>(&_header), sizeof(_header));
		_firstDataPos = _outputStream->tellp();
		_clustersInFile = 0;
		_writtenSourceNames.clear();
  }
  
  /**
	 * Open the file. After this call, ReadClusterMetaInfo() should be
	 * called to move the read position past the meta info. Once that has been
	 * done, the solution reading methods can be called.
	 */
	void OpenForReading(const char *filename)
	{
		delete _inputStream;
		_inputStream = new std::ifstream(filename);
		if(!_inputStream->good())
			throw std::runtime_error("Error reading input ionospheric solutions file");
		_inputStream->read(reinterpret_cast<char*>(&_header), sizeof(_header));
		if(_header.fileType != 2)
			throw std::runtime_error("Error reading ionospheric solution file; old format or file damaged");
		_clustersInFile = 0;
	}

  void ReadClusterMetaInfo(std::string& clusterName, std::vector<std::string>& sourceNames)
	{
		clusterName = readString();
		uint32_t count;
		_inputStream->read(reinterpret_cast<char*>(&count), sizeof(uint32_t));
		sourceNames.resize(count);
		for(std::vector<std::string>::iterator i=sourceNames.begin(); i!=sourceNames.end(); ++i)
			*i = readString();
		_firstDataPos = _inputStream->tellg();
		_clustersInFile++;
	}
	
	void ReadClusterMetaInfo(Model& expectedSources, std::vector<std::vector<ModelSource*>>& sourcesPerDirection)
	{
		_clustersInFile = 0;
		_inputStream->seekg(sizeof(_header), std::ios::beg);
		sourcesPerDirection.clear();
		std::map<std::string,size_t> sourceNameToIndex;
		for(size_t s=0; s!=expectedSources.SourceCount(); ++s)
		{
			if(sourceNameToIndex.count(expectedSources.Source(s).Name()))
				throw std::runtime_error("Model does not have unique source names");
			sourceNameToIndex.insert(std::make_pair(expectedSources.Source(s).Name(), s));
		}
		
		for(size_t d=0; d!=DirectionCount(); ++d)
		{
			std::string clusterName;
			std::vector<std::string> sourceNames;
			ReadClusterMetaInfo(clusterName, sourceNames);
			std::vector<ModelSource*> sources;
			for(size_t i=0; i!=sourceNames.size(); ++i)
			{
				std::map<std::string,size_t>::iterator sIter = sourceNameToIndex.find(sourceNames[i]);
				if(sIter == sourceNameToIndex.end()) throw std::runtime_error("Source not found in model");
				sources.push_back(&expectedSources.Source(sIter->second));
			}
			sourcesPerDirection.push_back(sources);
		}
	}

	double ReadSolution(IonSolutionType type, size_t interval, size_t channelBlock, size_t polarization, size_t direction)
	{
		Solution solution;
		ReadSolution(solution, interval, channelBlock, polarization, direction);
		if(!std::isfinite(solution.gain) || !std::isfinite(solution.dl) || !std::isfinite(solution.dm))
		{
			return std::numeric_limits<double>::quiet_NaN();
		}
		else {
			switch(type)
			{
				default:
				case GainSolution: return solution.gain;
				case DlSolution: return solution.dl;
				case DmSolution: return solution.dm;
			}
		}
	}
	
	double ReadAverageSolution(IonSolutionType type, size_t polarization, size_t direction)
	{
		double sum = 0.0;
		size_t count = 0;
		for(size_t i=0; i!=_header.intervalCount; ++i)
		{
			for(size_t c=0; c!=_header.channelBlockCount; ++c)
			{
				double solution = ReadSolution(type, i, c, polarization, direction);
				if(std::isfinite(solution))
				{
					sum += solution;
					++count;
				}
			}
		}
		return sum / count;
	}
	
  void ReadSolution(Solution& solution, size_t interval, size_t channelBlock, size_t polarization, size_t direction)
	{
		if(_clustersInFile != _header.directionCount)
			throw std::runtime_error("ReadSolution() called before all cluster meta data were read");
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelBlockCount + channelBlock) * _header.polarizationCount + polarization) * _header.directionCount + direction;
		_inputStream->seekg(_firstDataPos + sizeof(Solution) * index, std::ios::beg);
		if(!_inputStream->good())
			throw std::runtime_error("Error while reading solution from solutionfile: file corrupted?");
		_inputStream->read(reinterpret_cast<char*>(&solution), sizeof(Solution));
  }
  
  void WriteClusterMetaInfo(const std::string& clusterName, const std::vector<std::string>& sourceNames)
	{
		writeString(clusterName);
		uint32_t count = sourceNames.size();
		_outputStream->write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
		for(std::vector<std::string>::const_iterator i=sourceNames.begin(); i!=sourceNames.end(); ++i)
		{
			if(_writtenSourceNames.count(*i)>0)
				throw std::runtime_error("All sources in the model should have unique names!");
			_writtenSourceNames.insert(*i);
			writeString(*i);
		}
		_firstDataPos = _outputStream->tellp();
		++_clustersInFile;
	}

  void WriteSolution(const Solution& solution, size_t interval, size_t channel, size_t polarization, size_t direction)
  {
		if(_clustersInFile != _header.directionCount)
			throw std::runtime_error("WriteSolution() called before all cluster meta data were written");
		
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelBlockCount + channel) * _header.polarizationCount + polarization) * _header.directionCount + direction;
		_outputStream->seekp(_firstDataPos + sizeof(Solution) * index, std::ios::beg);
		_outputStream->write(reinterpret_cast<const char*>(&solution), sizeof(Solution));
  }

  void WriteChannelBlock(const Solution* solutions, size_t interval, size_t channel, size_t polarization)
  {
		std::unique_lock<std::mutex> lock(_mutex);
		size_t index = ((interval * _header.channelBlockCount + channel) * _header.polarizationCount + polarization) * _header.directionCount;
		_outputStream->seekp(_firstDataPos + sizeof(Solution) * index, std::ios::beg);
		_outputStream->write(reinterpret_cast<const char*>(solutions), sizeof(Solution) * _header.directionCount);
  }

 private:
  struct {
    char intro[8];
    uint32_t fileType;
    uint32_t structureType;
    uint32_t intervalCount, antennaCount, channelBlockCount, polarizationCount;
		uint32_t directionCount, parameterCount;
		double startTime, endTime;
		double startFrequency, endFrequency;
  } _header;
  std::ofstream *_outputStream;
  std::ifstream *_inputStream;
	std::mutex _mutex;
	size_t _firstDataPos, _clustersInFile;
	std::set<std::string> _writtenSourceNames;
	
	void writeString(const std::string& str)
	{
		uint32_t length = str.size();
		_outputStream->write(reinterpret_cast<const char*>(&length), sizeof(uint32_t));
		_outputStream->write(str.c_str(), str.size());
	}
	
	std::string readString()
	{
		uint32_t length;
		_inputStream->read(reinterpret_cast<char*>(&length), sizeof(uint32_t));
		if(!_inputStream->good())
			throw std::runtime_error("Error reading string length from solution file");
		ao::uvector<char> buffer(length+1);
		_inputStream->read(buffer.data(), length);
		if(!_inputStream->good())
			throw std::runtime_error("Error reading string from solution file");
		buffer[length] = 0;
		return std::string(buffer.data());
	}
};

#endif
