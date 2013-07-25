#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <fstream>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <stdint.h>
#include <cstdlib>

#include "banddata.h"

#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

using namespace casa;

void unwrap(long double **phases, const long double* const* weights, size_t channelCount, size_t polarizationCount)
{
	for(size_t p=0; p!=polarizationCount; ++p)
	{
		double phase = 0.0;
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(weights[ch][p] > 0.0)
			{
				if(std::isfinite(phases[ch][p]))
				{
					while(phases[ch][p] > phase + M_PI)
						phases[ch][p] -= 2.0 * M_PI;
					while(phases[ch][p] < phase - M_PI)
						phases[ch][p] += 2.0 * M_PI;
					phase = phases[ch][p];
				} else
					throw std::runtime_error("There are NaN's in the solution!");
			}
		}
	}
}

void subtractReference(long double** phases, const long double* const* refPhases, size_t channelCount, size_t polarizationCount)
{
	for(size_t ch=0; ch!=channelCount; ++ch)
	{
		phases[ch][3] -= refPhases[ch][3];
		phases[ch][0] -= refPhases[ch][0];
	}
}

void addXXandYY(long double** phases, size_t channelCount, size_t polarizationCount)
{
	for(size_t ch=0; ch!=channelCount; ++ch)
	{
		phases[ch][0] += phases[ch][3];
	}
}

void calculateOriginOffset(const long double* const* phases, const long double* const* weights, const BandData &bandData, size_t pol, size_t polarizationCount, long double &alpha, long double &beta)
{
	size_t channelCount = bandData.ChannelCount();
	long double avgY = 0.0, avgX = 0.0;
	long double weightSum = 0.0;
	for(size_t ch=0; ch!=channelCount; ++ch)
	{
		avgX += bandData.ChannelFrequency(ch) * weights[ch][pol];
		avgY += phases[ch][pol] * weights[ch][pol];
		weightSum += weights[ch][pol];
	}
	if(weightSum != 0.0)
	{
		avgX /= weightSum;
		avgY /= weightSum;
	}
	
	long double ssxx = 0.0, ssxy = 0.0;
	for(size_t ch=0;ch!=channelCount;ch++)
	{
		long double
			diffx = (long double) bandData.ChannelFrequency(ch) - avgX,
			diffy = phases[ch][pol] - avgY;
		ssxx += diffx*diffx * weights[ch][pol];
		ssxy += diffx*diffy * weights[ch][pol];
	}
	beta = (ssxx!=0.0) ? ssxy / ssxx : 0.0;
	alpha = avgY - beta * avgX;
}

static int polCharToIndex(char polarizationChar)
{
	switch(polarizationChar)
	{
		case 'X': case 'x':
		case 'R': case 'r':
		case 'I': case 'i':
			return 0;
		case 'Y': case 'y':
		case 'L': case 'l':
			return 1;
		default:
			throw std::runtime_error(std::string("Unknown pol char: ") + polarizationChar);
	}
}

/** Read the mapping between antennas and correlator inputs. */
void parseInputConfig(const char *filename, const std::vector<long double> &timeDelays)
{
	std::ifstream file(filename);
	std::ofstream outFile("instr_config_new.txt");
 
	std::string line;
	size_t nFlaggedInput = 0;
	size_t index = 0;
  while(file.good())
	{
		std::getline(file, line);
    if(!line.empty() && line[0]!='#')
		{
			std::istringstream str(line);
			
			size_t antennaIndex;
			std::string cableLen;
			unsigned dummy, inputFlag, isFlagged;
			char polChar;
			
			str >> dummy >> antennaIndex >> polChar >> cableLen;
			if(str.fail())
				throw std::runtime_error("Failed scanning instrument configuration file in line: " + line);
			str >> inputFlag;
			if(str.fail())
				isFlagged = false;
			else
				isFlagged = inputFlag;
			
			if(isFlagged) ++nFlaggedInput;
			
			// decode the string with the cable length. for a prefix of "EL_" this means the value is an electrical length
			// not a physical one, so no velocity factor should be applied.
			double cableLenDelta;
			if(cableLen.substr(0, 3) == "EL_")
				cableLenDelta = std::atof(&(cableLen.c_str()[3]));
			else
				throw std::runtime_error("Not in electrical length");

			unsigned char polarizationIndex = polCharToIndex(polChar);
			size_t inputIndex = index;
			
			const long double delay = timeDelays[antennaIndex];
			long double deltaDelay = delay * SPEED_OF_LIGHT;
			
			if(isFlagged)
				outFile << inputIndex << '\t' << antennaIndex << '\t' << ((polarizationIndex==0) ? 'X' : 'Y') << "\tEL_" << cableLenDelta << " 1 # skipped\n";
			else
				outFile << inputIndex << '\t' << antennaIndex << '\t' << ((polarizationIndex==0) ? 'X' : 'Y') << "\tEL_" << (cableLenDelta - deltaDelay) << " 0 # modified\n";
			
			++index;
		}
  }
  std::cout << "Read " << index << " inputs from " << filename << ", of which " << nFlaggedInput << " were flagged.\n";
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: fitdelays <ms> <phase-txt-file> <reference>\n"
			"-r will reverse the phases thus undo a previous run.\n";
	} else {
		size_t argi = 1;
		MeasurementSet ms(argv[argi]);
		std::ifstream phaseFile(argv[argi+1]);
		size_t referenceAntenna = atoi(argv[argi+2]);
		
		/**
		 * Read some meta data from the measurement set
		 */
		std::cout << "Opening measurement set..." << std::flush;
		MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		BandData bandData(spwTable);
		size_t channelCount = bandData.ChannelCount();
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ROArrayColumn<Complex> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ROArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		std::cout << " DONE\n";
		
		/**
		 * Read the phase offset file
		 */
		uint64_t antInFile, chanInFile, polInFile;
		phaseFile.read(reinterpret_cast<char*>(&antInFile), sizeof(antInFile));
		phaseFile.read(reinterpret_cast<char*>(&chanInFile), sizeof(chanInFile));
		phaseFile.read(reinterpret_cast<char*>(&polInFile), sizeof(polInFile));
		if(antInFile != antennaCount) throw std::runtime_error("Antenna counts do not match");
		if(polInFile != polarizationCount) throw std::runtime_error("Polarization counts do not match");
		if(channelCount%chanInFile!=0) throw std::runtime_error("Channel counts do not match");
		if(channelCount != chanInFile) throw std::runtime_error("These solutions were averaged: can't fit on averaged solutions");
		
		long double **phaseOffsets[antennaCount], **weights[antennaCount];
		if(chanInFile == 1)
		{
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				phaseOffsets[a] = new long double*[chanInFile];
				phaseOffsets[a][0] = new long double[polarizationCount];
				weights[a] = new long double*[chanInFile];
				weights[a][0] = new long double[polarizationCount];
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					double curPhaseOffset;
					phaseFile.read(reinterpret_cast<char*>(&curPhaseOffset), sizeof(curPhaseOffset));
					phaseOffsets[a][0][p] = curPhaseOffset;
					weights[a][0][p] = 0.0;
				}
			}
		} else {
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				phaseOffsets[a] = new long double*[chanInFile];
				weights[a] = new long double*[chanInFile];
				for(size_t ch=0; ch!=chanInFile; ++ch)
				{
					phaseOffsets[a][ch] = new long double[polarizationCount];
					weights[a][ch] = new long double[polarizationCount];
				}
			}
			for(size_t ch = 0; ch!=chanInFile; ++ch)
			{
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						double curPhaseOffset;
						phaseFile.read(reinterpret_cast<char*>(&curPhaseOffset), sizeof(curPhaseOffset));
						phaseOffsets[a][ch][p] = curPhaseOffset;
						weights[a][ch][p] = 0.0;
					}
				}
			}
		}
		
		std::cout << "Reading weights... " << std::flush;
		casa::Array<bool> flagArray(dataShape);
		for(size_t i=0; i!=ms.nrow(); ++i)
		{
			flagColumn.get(i, flagArray);
			size_t a1 = ant1Column(i), a2 = ant2Column(i);
			bool *flagArrayPtr = flagArray.cbegin();
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					if(!*flagArrayPtr)
					{
						weights[a1][ch][p]++;
						weights[a2][ch][p]++;
					}
					++flagArrayPtr;
				}
			}
		}
		
		std::vector<long double> weightsPerAntenna(antennaCount), weightsPerChannel(channelCount);
		for(size_t a=0; a!=antennaCount; ++a)
		{
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				weightsPerAntenna[a] += weights[a][ch][0];
				weightsPerAntenna[a] += weights[a][ch][3];
				weightsPerChannel[ch] += weights[a][ch][0];
				weightsPerChannel[ch] += weights[a][ch][3];
			}
		}
		std::cout << "DONE\nThe following channels have zero weight:";
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(fabs(weightsPerChannel[ch]) < 0.001)
				std::cout << ' ' << ch;
		}
		std::cout << "\nThe following antennas have zero weight:";
		for(size_t a=0; a!=antennaCount; ++a)
		{
			if(fabs(weightsPerAntenna[a]) < 0.001)
				std::cout << ' ' << a;
		}		
		std::cout << "\nSubtract reference antenna... " << std::flush;
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			if(a != referenceAntenna)
			{
				subtractReference(phaseOffsets[a], phaseOffsets[referenceAntenna], channelCount, polarizationCount);
			}
		}
		subtractReference(phaseOffsets[referenceAntenna], phaseOffsets[referenceAntenna], channelCount, polarizationCount);
		
		std::cout << "DONE\nAdding XX and YY... " << std::flush;
		for(size_t a = 0; a!=antennaCount; ++a)
			addXXandYY(phaseOffsets[a], channelCount, polarizationCount);
		
		std::cout << "DONE\nUnwrapping phases... " << std::flush;
		// Unwrap phases
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			unwrap(phaseOffsets[a], weights[a], channelCount, polarizationCount);
		}
		
		std::cout << "DONE\nFitting slopes... " << std::flush;
		// Find amount of phase wrapping towards origin
		std::vector<long double> timeDelays(antennaCount);
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			long double alpha, beta;
			calculateOriginOffset(phaseOffsets[a], weights[a], bandData, 0, polarizationCount, alpha, beta);
			
			// dphase = 2 pi nu dt
			// We found
			// dphase =      nu beta
			// So dt = beta / (2 pi)
			// However, remember this is XX and YY combined; need to divide by 2 once more.
			timeDelays[a] = beta / (4.0 * M_PI);
		}
		std::cout << "DONE\n";
		
		long double average = 0.0;
		for(size_t a = 0; a!=antennaCount; ++a)
			average += timeDelays[a];
		
		average /= antennaCount;
		std::cout << "Average time offset from reference antenna " << referenceAntenna << " : " << (1000000000.0 * average) << " ns\n";
		
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			timeDelays[a] -= average;
			std::cout << a << '\t' << timeDelays[a]*1000000000.0 << " ns\t" << timeDelays[a] * SPEED_OF_LIGHT << " m \n";
		}
		
		parseInputConfig("instr_config.txt", timeDelays);
	}
}
