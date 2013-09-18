#ifndef SOLUTION_APPLIER_H
#define SOLUTION_APPLIER_H

#include <complex>
#include <iostream>

#include <ms/MeasurementSets/MSAntenna.h>
#include <ms/MeasurementSets/MeasurementSet.h>

#include "banddata.h"
#include "solutionfile.h"
#include "matrix2x2.h"

class SolutionApplier
{
public:
	SolutionApplier() : _preset(false)
	{
	}
	
	void SetPresets(std::complex<double> xx, std::complex<double> xy, std::complex<double> yx, std::complex<double> yy)
	{
		_preset = true;
		_presetValues[0] = xx;
		_presetValues[1] = xy;
		_presetValues[2] = yx;
		_presetValues[3] = yy;
	}
	
	void Apply(casa::MeasurementSet& ms, SolutionFile& solutionFile)
	{
		/**
		 * Read some meta data from the measurement set
		 */
		std::cout << "Opening measurement set..." << std::flush;
		ms.reopenRW();
		casa::MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ArrayColumn<complex_t> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
		casa::ArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		if(polarizationCount != 4)
		  throw std::runtime_error("Should have 4 pols");
		
		std::cout << "DONE\nCounting timesteps... " << std::flush;
		double time = -1.0;
		std::vector<size_t> timestepRows;
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			if(timeColumn(rowIndex) != time)
			{
				timestepRows.push_back(rowIndex);
				time = timeColumn(rowIndex);
			}
		}
		size_t timestepCount = timestepRows.size();
		timestepRows.push_back(ms.nrow());
		std::cout << "DONE (" << timestepCount << " timesteps)\n";
	
		/**
		 * Read the solutions file
		 */
		std::vector<std::complex<double>*> values(antennaCount);
		for(size_t a = 0; a!=antennaCount; ++a) {
			values[a] = new std::complex<double>[channelCount*4];
		}
		if(_preset)
		{
			for(size_t a = 0; a!=antennaCount; ++a) {
				for(size_t ch = 0; ch!=channelCount; ++ch) {
					values[a][ch*4+0] = _presetValues[0];
					values[a][ch*4+1] = _presetValues[1];
					values[a][ch*4+2] = _presetValues[2];
					values[a][ch*4+3] = _presetValues[3];
				}		  
			}
		}
		else {
			std::cout << "Checking solutions file..." << std::flush;
			if(solutionFile.AntennaCount() != antennaCount) throw std::runtime_error("Antenna counts do not match");
			if(solutionFile.PolarizationCount() != polarizationCount) throw std::runtime_error("Polarization counts do not match");
			if(channelCount%solutionFile.ChannelCount()!=0) throw std::runtime_error("Channel counts do not match");
			std::cout << " DONE\n";
		}
		
		/**
		 * Apply corrections
		 */
		std::cout << "Applying solutions...\n";
		casa::Array<complex_t> data(dataShape);
		for(size_t interval=0; interval!=solutionFile.IntervalCount(); ++interval)
		{
			// Read the solutions for this interval
			if(!_preset)
			{
				for(size_t a = 0; a!=antennaCount; ++a) {
					for(size_t ch = 0; ch!=channelCount; ++ch) {
						for(size_t p = 0; p!=4; ++p) {
							values[a][ch*4+p] = solutionFile.ReadNextSolution();
						}
					}
				}
			}
			
			size_t
				intervalTimestepStart = (interval*timestepCount) / solutionFile.IntervalCount(),
				intervalTimestepEnd = ((interval+1)*timestepCount) / solutionFile.IntervalCount(),
				intervalRowStart = timestepRows[intervalTimestepStart],
				intervalRowEnd = timestepRows[intervalTimestepEnd];
			std::cout << "- Interval " << (interval+1) << '/' << solutionFile.IntervalCount() << " (" << intervalRowStart << '-' << intervalRowEnd << ")\n";
			std::cout << "  Antenna1: " << values[1][72*4] << "\n";
			for(size_t rowIndex=intervalRowStart; rowIndex!=intervalRowEnd; ++rowIndex)
			{
				// Cross correlation?
				size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
				if(a1 != a2) {
					dataColumn.get(rowIndex, data);
					casa::Array<complex_t>::contiter dataPtr = data.cbegin();
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						size_t chFileIndex = ch * 4;
						std::complex<double>
						*solA = &values[a1][chFileIndex],
						*solB = &values[a2][chFileIndex];
						std::complex<double> dataVals[4] = {
							dataPtr[0], dataPtr[1], dataPtr[2], dataPtr[3]
						};
						applySolution(dataVals, solA, solB);
						dataPtr[0] = dataVals[0];
						dataPtr[1] = dataVals[1];
						dataPtr[2] = dataVals[2];
						dataPtr[3] = dataVals[3];
						dataPtr += 4;
					}
				}
				dataColumn.put(rowIndex, data);
			}
		}
		
		/**
		 * Free mem
		 */
		for(size_t a = 0; a!=antennaCount; ++a)
		{
		  delete[] values[a];
		}
	}
private:
	void applySolution(std::complex<double> *dataVal, const std::complex<double> *solA, const std::complex<double> *solB)
	{
		std::complex<double> solATimesData[4];
		Matrix2x2::ATimesB(solATimesData, solA, dataVal);
		Matrix2x2::ATimesHermB(dataVal, solATimesData, solB);
	}

	bool _preset;
	std::complex<double> _presetValues[4];
};

#endif
