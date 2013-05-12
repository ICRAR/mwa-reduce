#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include "stdint.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

using namespace casa;

void PhaseRotate(std::complex<float> &dest, long double phaseOffset)
{
	long double phCos = cosl(phaseOffset), phSin = sinl(phaseOffset);
	long double r = dest.real(),	i = dest.imag();
	dest = std::complex<float>(
		r * phCos - i * phSin,
		r * phSin + i * phCos
	);
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: applyphaseoffset [-r] <ms> <phase-txt-file>\n"
			"-r will reverse the phases thus undo a previous run.\n";
	} else {
		bool reverse = false;
		size_t argi = 1;
		if(strcmp(argv[argi], "-r")==0)
		{
			reverse=true;
			++argi;
		}
		MeasurementSet ms(argv[argi], Table::Update);
		std::ifstream phaseFile(argv[argi+1]);
		
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
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		std::cout << " DONE\n";
		
		/**
		 * Read the phase offset file
		 */
		std::cout << "Reading phase offset file..." << std::flush;
		uint64_t antInFile, chanInFile, polInFile;
		phaseFile.read(reinterpret_cast<char*>(&antInFile), sizeof(antInFile));
		phaseFile.read(reinterpret_cast<char*>(&chanInFile), sizeof(chanInFile));
		phaseFile.read(reinterpret_cast<char*>(&polInFile), sizeof(polInFile));
		if(antInFile != antennaCount) throw std::runtime_error("Antenna counts do not match");
		if(polInFile != polarizationCount) throw std::runtime_error("Polarization counts do not match");
		if(channelCount%chanInFile!=0) throw std::runtime_error("Channel counts do not match");
		
		long double **phaseOffsets[antennaCount];
		if(chanInFile == 1)
		{
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				phaseOffsets[a] = new long double*[chanInFile];
				phaseOffsets[a][0] = new long double[polarizationCount];
				//size_t antNo;
				//phaseFile >> antNo;
				//if(antNo != a) throw std::runtime_error("File corrupted: antenna numbers do not match");
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					double curPhaseOffset;
					phaseFile.read(reinterpret_cast<char*>(&curPhaseOffset), sizeof(curPhaseOffset));
					phaseOffsets[a][0][p] = curPhaseOffset;
				}
			}
		} else {
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				phaseOffsets[a] = new long double*[chanInFile];
				for(size_t ch=0; ch!=chanInFile; ++ch)
					phaseOffsets[a][ch] = new long double[polarizationCount];
			}
			for(size_t ch = 0; ch!=chanInFile; ++ch)
			{
				//size_t chNo;
				//phaseFile >> chNo;
				//if(chNo != ch) throw std::runtime_error("File correupted: Channel numbers do no match");
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						double curPhaseOffset;
						phaseFile.read(reinterpret_cast<char*>(&curPhaseOffset), sizeof(curPhaseOffset));
						phaseOffsets[a][ch][p] = curPhaseOffset;
					}
				}
			}
		}
		
		if(reverse)
		{
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				for(size_t ch = 0; ch!=chanInFile; ++ch)
				{
					for(size_t p = 0; p!=polarizationCount; ++p)
						phaseOffsets[a][ch][p] = -phaseOffsets[a][ch][p];
				}
			}
		}
		
		std::cout << " DONE\n";
		
		/**
		 * Apply corrections
		 */
		std::cout << "Applying phase offsets..." << std::flush;
		Array<complex_t> data(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				long double sign;
				if(a1 > a2) {
					std::swap(a1, a2);
					sign = -1.0;
				} else {
					sign = 1.0;
				}
				dataColumn.get(rowIndex, data);
				Array<complex_t>::iterator dataPtr = data.begin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					size_t fileChIndex = ch*chanInFile/channelCount;
					long double *chanPhaseOffsets1, *chanPhaseOffsets2;
					chanPhaseOffsets1 = phaseOffsets[a1][fileChIndex];
					chanPhaseOffsets2 = phaseOffsets[a2][fileChIndex];
					if(polarizationCount == 1)
					{
						long double ph = chanPhaseOffsets2[0] - chanPhaseOffsets1[0];
						PhaseRotate(*dataPtr, ph*sign);
						++dataPtr;
					} else if(polarizationCount == 2)
					{
						long double
							phXX = chanPhaseOffsets2[0] - chanPhaseOffsets1[0],
							phYY = chanPhaseOffsets2[1] - chanPhaseOffsets1[1];
						PhaseRotate(*dataPtr, phXX*sign);
						++dataPtr;
						PhaseRotate(*dataPtr, phYY*sign);
						++dataPtr;
					} else if(polarizationCount == 4)
					{
						long double
							phXX = chanPhaseOffsets2[0] - chanPhaseOffsets1[0],
							phXY = chanPhaseOffsets2[3] - chanPhaseOffsets1[0],
							phYX = chanPhaseOffsets2[0] - chanPhaseOffsets1[3],
							phYY = chanPhaseOffsets2[3] - chanPhaseOffsets1[3];
						PhaseRotate(*dataPtr, phXX*sign);
						++dataPtr;
						PhaseRotate(*dataPtr, phXY*sign);
						++dataPtr;
						PhaseRotate(*dataPtr, phYX*sign);
						++dataPtr;
						PhaseRotate(*dataPtr, phYY*sign);
						++dataPtr;
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
			for(size_t ch=0; ch!=chanInFile; ++ch)
				delete[] phaseOffsets[a][ch];
			delete[] phaseOffsets[a];
		}
			
		std::cout << " DONE\n";		
	}
}
