#include "solutionfile.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "stdint.h"

using namespace casa;

void ApplySolution(std::complex<double> *dataVal, const std::complex<double> *solA, const std::complex<double> *solB)
{
  std::complex<double> solATimesData[4];
  solATimesData[0] = solA[0] * dataVal[0] + solA[1] * dataVal[2];
  solATimesData[1] = solA[0] * dataVal[1] + solA[1] * dataVal[3];
  solATimesData[2] = solA[2] * dataVal[0] + solA[3] * dataVal[2];
  solATimesData[3] = solA[2] * dataVal[1] + solA[3] * dataVal[3];

  dataVal[0] = solATimesData[0] * std::conj(solB[0]) + solATimesData[1] * std::conj(solB[1]);
  dataVal[1] = solATimesData[0] * std::conj(solB[2]) + solATimesData[1] * std::conj(solB[3]);
  dataVal[2] = solATimesData[2] * std::conj(solB[0]) + solATimesData[3] * std::conj(solB[1]);
  dataVal[3] = solATimesData[2] * std::conj(solB[2]) + solATimesData[3] * std::conj(solB[3]);
}

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std::cout << "Usage: applysolutions [s xx xy yx yy] <ms> <gains-bin-file>\n"
				"Will apply the found solution matrices.\n";
    } else {
		size_t argi = 1;
		double xx=0.0, xy=0.0, yx=0.0, yy=0.0;
		bool preset = false;
		if(strcmp(argv[argi], "-s") == 0)
		{
			xx = atof(argv[argi+1]);
			xy = atof(argv[argi+2]);
			yx = atof(argv[argi+3]);
			yy = atof(argv[argi+4]);
			argi += 5;
			preset = true;
		}
		MeasurementSet ms(argv[argi], Table::Update);
		/**
		 * Read some meta data from the measurement set
		 */
		std::cout << "Opening measurement set..." << std::flush;
		MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		if(polarizationCount != 4)
		  throw std::runtime_error("Should have 4 pols");
		std::cout << " DONE\n";
		
		/**
		 * Read the solutions file
		 */
		std::vector<std::complex<double>*> values(antennaCount);
		for(size_t a = 0; a!=antennaCount; ++a) {
			values[a] = new std::complex<double>[channelCount*4];
		}
		if(preset)
		{
			for(size_t a = 0; a!=antennaCount; ++a) {
				for(size_t ch = 0; ch!=channelCount; ++ch) {
					values[a][ch*4+0] = xx;
					values[a][ch*4+1] = xy;
					values[a][ch*4+2] = yx;
					values[a][ch*4+3] = yy;
				}		  
			}
		}
		else {
			std::cout << "Reading solutions file..." << std::flush;
			SolutionFile solutionFile;
			solutionFile.OpenForReading(argv[2]);
			if(solutionFile.AntennaCount() != antennaCount) throw std::runtime_error("Antenna counts do not match");
			if(solutionFile.PolarizationCount() != polarizationCount) throw std::runtime_error("Polarization counts do not match");
			if(channelCount%solutionFile.ChannelCount()!=0) throw std::runtime_error("Channel counts do not match");
			
			for(size_t a = 0; a!=antennaCount; ++a) {
				for(size_t ch = 0; ch!=channelCount; ++ch) {
					for(size_t p = 0; p!=4; ++p) {
						values[a][ch*4+p] = solutionFile.ReadNextSolution();
					}
				}		  
			}
			std::cout << " DONE\n";
		}
		
		/**
		 * Apply corrections
		 */
		std::cout << "Applying solutions..." << std::flush;
		Array<complex_t> data(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2) {
				dataColumn.get(rowIndex, data);
				Array<complex_t>::contiter dataPtr = data.cbegin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					size_t chFileIndex = ch * 4;
					std::complex<double>
					*solA = &values[a1][chFileIndex],
					*solB = &values[a2][chFileIndex];
					std::complex<double> dataVals[4] = {
						dataPtr[0], dataPtr[1], dataPtr[2], dataPtr[3]
					};
					ApplySolution(dataVals, solA, solB);
					dataPtr[0] = dataVals[0];
					dataPtr[1] = dataVals[1];
					dataPtr[2] = dataVals[2];
					dataPtr[3] = dataVals[3];
					dataPtr += 4;
				}
			}
		  dataColumn.put(rowIndex, data);
		}
		
		/**
		 * Free mem
		 */
		for(size_t a = 0; a!=antennaCount; ++a)
		{
		  delete[] values[a];
		}
		
		std::cout << " DONE\n";
	}
}
