#include <iostream>

#include <string.h>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <stdexcept>

#include <cmath>
#include <fstream>

#include "sourcesdf.h"
#include "banddata.h"

using namespace casa;

typedef std::complex<float> complex_t;

std::string GainToString(long double gain)
{
	std::stringstream s;
	const char *unit;
	if(gain >= 1.0)
	{
		unit = "Jy";
	} else if(gain >= 0.001)
	{
		gain *= 1000.0;
		unit = "mJy";
	} else
	{
		gain *= 1000000.0;
		unit = "μJy";
	}
	s.precision(3);
	s << gain << ' ' << unit;
	return s.str();
}

size_t ValueIndex(size_t a1, size_t a2, size_t antennaCount)
{
	return (a1>a2) ? (a2*antennaCount + a1) : (a1*antennaCount + a2);
}
size_t ValueIndexA2highest(size_t a1, size_t a2, size_t antennaCount)
{
	return a1*antennaCount + a2;
}
size_t ValueIndexA1highest(size_t a1, size_t a2, size_t antennaCount)
{
	return a2*antennaCount + a1;
}

void DivideVisibilitiesByCount(long double **visValues, const size_t* const* counts, size_t antennaCount, size_t polarizationCount, size_t channelCount)
{
	for(size_t a1=0;a1!=antennaCount;++a1)
	{
		for(size_t a2=a1+1;a2!=antennaCount;++a2)
		{
			size_t index = ValueIndexA2highest(a1, a2, antennaCount);
			long double *antValues = visValues[index];
			const size_t *antCounts = counts[index];
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					if(*antCounts == 0)
						*antValues = 0.0;
					else
						*antValues = (*antValues) / *antCounts;
					++antValues;
					++antCounts;
				}
			}
		}
	}
}

void InitializeWeights(long double **visWeights, const size_t* const* counts, size_t antennaCount, size_t polarizationCount, size_t channelCount, long double maxCount)
{
	const long double normFactor = 1.0 / maxCount;
	for(size_t a1=0;a1!=antennaCount;++a1)
	{
		for(size_t a2=a1+1;a2!=antennaCount;++a2)
		{
			size_t index = ValueIndexA2highest(a1, a2, antennaCount);
			long double *antWeights = visWeights[index];
			const size_t *antCounts = counts[index];
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					*antWeights = (long double) (*antCounts) * normFactor;
					++antWeights;
					++antCounts;
				}
			}
		}
	}
}

size_t MaximumCount(const size_t* const* counts, size_t antennaCount, size_t polarizationCount, size_t channelCount)
{
	size_t maxCount = 0;
	for(size_t a1=0;a1!=antennaCount;++a1)
	{
		for(size_t a2=a1+1;a2!=antennaCount;++a2)
		{
			size_t index = ValueIndexA2highest(a1, a2, antennaCount);
			const size_t *antCounts = counts[index];
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					if(*antCounts > maxCount) maxCount = *antCounts;
					++antCounts;
				}
			}
		}
	}
	return maxCount;
}

long double AntennaGain(size_t antenna, size_t antennaCount, size_t eIndex, const long double* const* visValues, long double refVis, const long double* const* visWeights)
{
	long double gain = 0.0, weight = 0.0;
	for(size_t a2 = 0; a2 != antenna; ++a2)
	{
		size_t aIndex = ValueIndex(antenna, a2, antennaCount);
		long double curWeight = visWeights[aIndex][eIndex];
		gain += visValues[aIndex][eIndex] * curWeight;
		weight += curWeight;
	}
	for(size_t a2 = antenna+1; a2 != antennaCount; ++a2)
	{
		size_t aIndex = ValueIndex(antenna, a2, antennaCount);
		long double curWeight = visWeights[aIndex][eIndex];
		gain += visValues[aIndex][eIndex] * curWeight;
		weight += curWeight;
	}
	if(weight != 0.0)
		gain /= (weight * refVis);
	else
		gain = 1.0;
	return gain;
}

long double AntennaWeight(size_t antenna, size_t antennaCount, size_t eIndex, const long double* const* visWeights)
{
	long double weight = 0;
	for(size_t a2 = 0; a2 != antenna; ++a2)
	{
		size_t aIndex = ValueIndex(antenna, a2, antennaCount);
		weight += visWeights[aIndex][eIndex];
	}
	for(size_t a2 = antenna+1; a2 != antennaCount; ++a2)
	{
		size_t aIndex = ValueIndex(antenna, a2, antennaCount);
		weight += visWeights[aIndex][eIndex];
	}
	weight /= (long double) (antennaCount-1);
	return weight;
}

std::pair<long double, long double> ChannelStdError(size_t ch, size_t antennaCount, size_t polarizationCount, const long double* const* visValues, long  double sourceFlux, const long double* const* visWeights)
{
	long double stdError = 0.0, weightSum = 0.0;
	size_t eIndex = polarizationCount * ch;
	for(size_t pol=0; pol!=polarizationCount; ++pol)
	{
		if(polarizationCount!=4 || pol==0 || pol==3) {
			for(size_t a1=0;a1!=antennaCount;++a1)
			{
				long double gain = AntennaGain(a1, antennaCount, eIndex, visValues, sourceFlux, visWeights);
				
				if(polarizationCount!=4 || pol==0 || pol==3)
				{
					long double weight = AntennaWeight(a1, antennaCount, eIndex, visWeights);
					stdError += (gain-1.0) * (gain-1.0) * weight;
					weightSum += weight;
				}
			}
		}
	}
	if(weightSum == 0.0)
		return std::pair<long double, long double>(0.0, 0.0);
	else
		return std::pair<long double, long double>(sqrtl(stdError / weightSum), weightSum);
}

int main(int argc, char *argv[])
{
	if(argc < 6)
	{
		std::cout
			<< "Usage: calcgains [-a] <measurementset.ms> <source-flux> <spectral index> <ref frequency> <output.txt>\n\n"
			<< "This will calculate \"static\" gains for all stations. To do this, it assumes the array\n"
			<< "is tracking a bright source in the phase centre. It produces approximate least-squares solutions.\n"
			<< "Option -a will average over frequency before fitting.\n"
			<< "Frequency given in MHz, flux density in Jy.\n";
	} else {
		int argi = 1;
		size_t avgChannelCount = 0;
		while(argi<argc && argv[argi][0]=='-')
		{
			if(strcmp(argv[argi], "-a")==0) {
				++argi;
				avgChannelCount = atoi(argv[argi]);
				++argi;
			}
		}
		const char *msName = argv[argi];
		const long double
			sourceFluxDensity = atof(argv[argi+1]),
			spectralIndex = atof(argv[argi+2]),
			refFrequency = atof(argv[argi+3]);
		const char *outName = argv[argi+4];
		const SourceSDFWithSI<long double>
			sourceSDF(sourceFluxDensity, spectralIndex, refFrequency*1000000.0);
		
		MeasurementSet ms(msName);
		
		std::cout << "Reading measurement set... " << std::flush;
		
		/**
		 * Read some meta data from the measurement set
		 */
		MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		if(avgChannelCount == 0) avgChannelCount = channelCount;
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ROArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		/**
		 * Initialize memory
		 */
		const size_t matrixSize = antennaCount * antennaCount;
		size_t sumsPerBaseline = channelCount * polarizationCount;
		long double *visValues[matrixSize], *visWeights[matrixSize];
		size_t *counts[matrixSize];
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			visValues[e] = new long double[sumsPerBaseline];
			visWeights[e] = new long double[sumsPerBaseline];
			for(size_t i=0; i<sumsPerBaseline; ++i) visValues[e][i] = 0.0;
			counts[e] = new size_t[sumsPerBaseline];
			memset(&counts[e][0], 0, sizeof(size_t) * sumsPerBaseline);
		}
		
		long double *gains[antennaCount];
		for(size_t a=0; a!=antennaCount; ++a) 
		{
			gains[a] = new long double[sumsPerBaseline];
			for(size_t i=0; i!=sumsPerBaseline; ++i) gains[a][i] = 1.0;
		}
		
		/**
		 * Read the average real parts of the visibilities per baseline per channel.
		 * Average channels if requested
		 */
		{
			Array<complex_t> data(dataShape);
			Array<bool> flags(dataShape);
			for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
			{
				// Cross correlation?
				size_t antenna1 = ant1Column.get(rowIndex), antenna2 = ant2Column.get(rowIndex);
				if(antenna1 != antenna2)
				{
					dataColumn.get(rowIndex, data);
					flagColumn.get(rowIndex, flags);
					Array<complex_t>::const_iterator dataPtr = data.begin();
					Array<bool>::const_iterator flagPtr = flags.begin();
					size_t index = ValueIndex(antenna1, antenna2, antennaCount);
					long double *baselineValues = visValues[index];
					size_t *baselineCounts = counts[index];
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						size_t eIndex = (ch*avgChannelCount/channelCount) * polarizationCount;
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							// Only include non-flagged data
							if(!(*flagPtr) && std::isfinite(dataPtr->real()) && std::isfinite(dataPtr->imag())) {
								baselineValues[eIndex] += dataPtr->real();
								baselineCounts[eIndex]++;
							}
							++dataPtr;
							++flagPtr;
							++eIndex;
						}
					}
				}
			}
		}
		
		std::cout << "DONE!\n";
		
		/**
		 * Apply spectral index to source flux
		 */
		double sourceFlux[avgChannelCount];
		for(size_t ch=0;ch!=avgChannelCount;++ch)
		{
			size_t beginIndex = ch*channelCount/avgChannelCount;
			size_t endIndex = (ch+1)*channelCount/avgChannelCount-1;
			if(beginIndex == endIndex)
				sourceFlux[ch] = sourceSDF.FluxAtFrequency(bandData.ChannelFrequency(beginIndex));
			else
			{
				double beginFreq = bandData.ChannelFrequency(beginIndex);
				double endFreq = bandData.ChannelFrequency(endIndex);
				sourceFlux[ch] = sourceSDF.IntegratedFlux(beginFreq, endFreq);
			}
		}
		std::cout << "Source flux: " << GainToString(sourceFlux[0]) << " @ " <<
			bandData.ChannelFrequency(0)/1000000.0 << " MHz - "
			<< GainToString(sourceFlux[avgChannelCount-1]) << " @ " <<
			bandData.ChannelFrequency(channelCount-1)/1000000.0 << " MHz " << '\n';
		
		size_t maxCount = MaximumCount(counts, antennaCount, polarizationCount, channelCount);
		InitializeWeights(visWeights, counts, antennaCount, polarizationCount, avgChannelCount, maxCount);
		
		channelCount = avgChannelCount;
		sumsPerBaseline = polarizationCount * avgChannelCount;

		DivideVisibilitiesByCount(visValues, counts, antennaCount, polarizationCount, channelCount);
		
		/**
		 * Calculate initial error
		 */
		std::cout << "Initial standard error of gains: " << std::flush;
		long double stdError = 0.0, stdErrorWeight = 0.0;
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, polarizationCount, visValues, sourceFlux[ch], visWeights);
			stdError += result.first;
			stdErrorWeight += result.second;
		}
		stdError /= stdErrorWeight;
		std::cout << GainToString(stdError) << ".\n";
		
		size_t iteration=0, gotWorseCount=0;
		long double prevAbsError = 10.0, minError = 10.0, stepSize = 0.25;
		do
		{
			std::cout << "Iteration " << iteration << ": " << std::flush;
			prevAbsError = stdError;
			stdError = 0.0;
			stdErrorWeight = 0.0;
			
			/**
			* Calculate per antenna gains
			*/
			size_t eIndex = 0;
			for(size_t ch = 0; ch!=channelCount; ++ch)
			{
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					// Skip XY / YX if they exist
					if(polarizationCount!=4 || p==0 || p==3)
					{
						for(size_t a1 = 0; a1!=antennaCount; ++a1)
						{
							long double antennaGain = AntennaGain(a1, antennaCount, eIndex, visValues, sourceFlux[ch], visWeights);
							long double antennaWeight = AntennaWeight(a1, antennaCount, eIndex, visWeights);
							
							if(eIndex/polarizationCount == (channelCount/2) && a1==0)
							{
								std::cout << GainToString(antennaGain) << '\t' << GainToString(gains[a1][eIndex]) << '\t' << "W=" << antennaWeight << '\t';
							}
							
							// Apply the average gain to the 'gains' matrix (i.e., step towards desired gain)
							long double step = 1.0 + (antennaGain-1.0) * stepSize * antennaWeight;
							gains[a1][eIndex] *= step;
							if(step != 0.0)
							{
								for(size_t a2 = 0; a2!=a1; ++a2)
								{
									size_t aIndex = ValueIndexA1highest(a1, a2, antennaCount);
									visValues[aIndex][eIndex] /= step;
								}
								for(size_t a2 = a1+1; a2!=antennaCount; ++a2)
								{
									size_t aIndex = ValueIndexA2highest(a1, a2, antennaCount);
									visValues[aIndex][eIndex] /= step;
								}
							} else std::cout << "Ant " << a1 << " made step to 0.\n";
							
						}
					}
					++eIndex;
				}
				
				std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, polarizationCount, visValues, sourceFlux[ch], visWeights);
				stdError += result.first;
				stdErrorWeight += result.second;
			}
			stdError /= stdErrorWeight;
			if(stdError < minError) minError = stdError;
			std::cout << "Standard error: " << GainToString(stdError) << '\n';
			
			bool gotBetter = stdError < prevAbsError;
			if(!gotBetter) { ++gotWorseCount; stepSize*=0.9; }
			++iteration;
		} while((gotWorseCount < 10 || stdError > minError*1.1) && iteration < 100);
		
		/**
		* Print results
		*/
		std::ofstream outFile(outName);
		outFile.precision(10);
		outFile << antennaCount << '\t' << channelCount << '\t' << polarizationCount << '\n';
		if(channelCount == 1)
		{
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				outFile << a;
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					long double *antennaGains = gains[a];
					outFile << '\t' << antennaGains[p];
				}
				outFile << '\n';
			}
		} else {
			size_t eIndex = 0;
			for(size_t ch = 0; ch!=channelCount; ++ch)
			{
				outFile << ch;
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						long double *antennaGains = gains[a];
						outFile << '\t' << antennaGains[eIndex];
					}
					++eIndex;
				}
				outFile << '\n';
			}
		}
		
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			delete[] visValues[e];
			delete[] counts[e];
		}
		for(size_t a=0; a!=antennaCount;++a)
			delete[] gains[a];
	}
	
	return 0;
}
