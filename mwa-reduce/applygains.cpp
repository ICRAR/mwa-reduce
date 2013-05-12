#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "stdint.h"

using namespace casa;

void ApplyGain(std::complex<float> &val, bool &flag, long double g)
{
	if(g == 0.0)
		flag = true;
	else
		val = std::complex<float>(val.real()/g, val.imag()/g);
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: applygains <ms> <gains-bin-file>\n"
			"Will apply the found gains. Visibilities with zero gain will be flagged.";
	} else {
		MeasurementSet ms(argv[1], Table::Update);
		std::ifstream gainsFile(argv[2]);
		
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
		std::cout << " DONE\n";
		
		/**
		 * Read the gains file
		 */
		std::cout << "Reading gains file..." << std::flush;
		uint64_t antInFile, chanInFile, polInFile;
		gainsFile.read(reinterpret_cast<char*>(&antInFile), sizeof(antInFile));
		gainsFile.read(reinterpret_cast<char*>(&chanInFile), sizeof(chanInFile));
		gainsFile.read(reinterpret_cast<char*>(&polInFile), sizeof(polInFile));
		if(antInFile != antennaCount) throw std::runtime_error("Antenna counts do not match");
		if(polInFile != polarizationCount) throw std::runtime_error("Polarization counts do not match");
		if(channelCount%chanInFile!=0) throw std::runtime_error("Channel counts do not match");
		
		long double **gains[antennaCount];
		if(chanInFile == 1)
		{
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				gains[a] = new long double*[chanInFile];
				gains[a][0] = new long double[polarizationCount];
				//size_t antNo;
				//gainsFile >> antNo;
				//if(antNo != a) throw std::runtime_error("File corrupted: antenna numbers do not match");
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					double curGain;
					gainsFile.read(reinterpret_cast<char*>(&curGain), sizeof(curGain));
					gains[a][0][p] = curGain;
				}
			}
		} else {
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				gains[a] = new long double*[chanInFile];
				for(size_t ch=0; ch!=chanInFile; ++ch)
					gains[a][ch] = new long double[polarizationCount];
			}
			for(size_t ch = 0; ch!=chanInFile; ++ch)
			{
				//size_t chNo;
				//gainsFile >> chNo;
				//if(chNo != ch) throw std::runtime_error("File corrupted: Channel numbers do no match");
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						double curGain;
						gainsFile.read(reinterpret_cast<char*>(&curGain), sizeof(curGain));
						gains[a][ch][p] = curGain;
					}
				}
			}
		}
		std::cout << " DONE\n";
		
		/**
		 * Apply corrections
		 */
		std::cout << "Applying gains..." << std::flush;
		Array<complex_t> data(dataShape);
		Array<bool> flags(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				if(a1 > a2) std::swap(a1, a2);
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				Array<complex_t>::iterator dataPtr = data.begin();
				Array<bool>::iterator flagPtr = flags.begin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					long double *gains1, *gains2;
					size_t chFileIndex = ch * chanInFile / channelCount;
					gains1 = gains[a1][chFileIndex];
					gains2 = gains[a2][chFileIndex];
					if(polarizationCount == 1)
					{
						long double g = gains2[0] * gains1[0];
						ApplyGain(*dataPtr, *flagPtr, g);
						++dataPtr; ++flagPtr;
					} else if(polarizationCount == 2)
					{
						long double
							gXX = gains2[0] * gains1[0],
							gYY = gains2[1] * gains1[1];
						ApplyGain(*dataPtr, *flagPtr, gXX);
						++dataPtr; ++flagPtr;
						ApplyGain(*dataPtr, *flagPtr, gYY);
						++dataPtr; ++flagPtr;
					} else if(polarizationCount == 4)
					{
						long double
							gXX = gains2[0] * gains1[0],
							gXY = gains2[3] * gains1[0],
							gYX = gains2[0] * gains1[3],
							gYY = gains2[3] * gains1[3];
						ApplyGain(*dataPtr, *flagPtr, gXX);
						++dataPtr; ++flagPtr;
						ApplyGain(*dataPtr, *flagPtr, gXY);
						++dataPtr; ++flagPtr;
						ApplyGain(*dataPtr, *flagPtr, gYX);
						++dataPtr; ++flagPtr;
						ApplyGain(*dataPtr, *flagPtr, gYY);
						++dataPtr; ++flagPtr;
					}
				}
				dataColumn.put(rowIndex, data);
				flagColumn.put(rowIndex, flags);
			}
		}
		
		/**
		 * Free mem
		 */
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			for(size_t ch=0; ch!=chanInFile; ++ch)
				delete[] gains[a][ch];
			delete[] gains[a];
		}
		
		std::cout << " DONE\n";
	}
}
