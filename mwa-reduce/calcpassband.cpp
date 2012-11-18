#include <iostream>

#include <string.h>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <stdexcept>

#include <cmath>

using namespace casa;

template<typename T>
T stddev(T sum, T sumSq, size_t n)
{
	if(n == 0)
		return T(0.0);
	else {
		T sumMeanSquared = sum * sum / n;
		return sqrtl((sumSq - sumMeanSquared) / n);
	}
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: calcpassband [-n] <ms> <nrsubbands> <list of selected subband indices>\n-n will normalize the standard deviations so that the median is one.\n";
	} else {
		size_t argi = 1;
		bool doNormalize = false;
		if(strcmp(argv[argi], "-n") == 0)
		{
			++argi;
			doNormalize = true;
		}
		MeasurementSet ms(argv[argi]);
		
		size_t collapseCount = atoi(argv[argi+1]);
		
		std::vector<size_t> selectedSubBands;
		for(int i=argi+2;i!=argc;++i)
			selectedSubBands.push_back(atoi(argv[i]));
		
		//MSAntenna aTable = ms.antenna();
		//size_t antennaCount = aTable.nrow();
		
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if((channelCount%collapseCount) != 0) throw std::runtime_error("Number of channels is not divisable by given subband count");
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		size_t channelsPerSubband = channelCount / collapseCount;
		
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ROArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		size_t sumCount = channelsPerSubband * polarizationCount;
		std::complex<long double> sums[sumCount], sumsSq[sumCount];
		size_t counts[sumCount];
		memset(counts, 0, sizeof(size_t) * sumCount);
		
		Array<complex_t> data(dataShape);
		Array<bool> flags(dataShape);
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			// Cross correlation?
			if(ant1Column.get(rowIndex) != ant2Column.get(rowIndex))
			{
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				Array<complex_t>::const_iterator dataPtr = data.begin();
				Array<bool>::const_iterator flagPtr = flags.begin();
				for(size_t sb=0;sb<collapseCount;++sb)
				{
					// Is this sub-band not selected?
					if(std::find(selectedSubBands.begin(), selectedSubBands.end(), sb) == selectedSubBands.end())
					{
						for(size_t i=0;i<sumCount;++i) 	{ ++dataPtr; ++flagPtr; }
					}
					else {
						size_t sumIndex = 0;
						for(size_t i=0;i!=channelsPerSubband;++i)
						{
							for(size_t polarizationIndex=0;polarizationIndex!=polarizationCount;++polarizationIndex)
							{
								if(!(*flagPtr) && std::isfinite(dataPtr->real()) && std::isfinite(dataPtr->imag()))
								{
									sums[sumIndex] += *dataPtr;
									sumsSq[sumIndex] += std::complex<long double>(dataPtr->real() * dataPtr->real(), dataPtr->imag() * dataPtr->imag());
									++counts[sumIndex];
								}
								++dataPtr;
								++flagPtr;
								++sumIndex;
							}
						}
					}
				}
			}
		}
		
		// We will make the passband median "1", so that we do not scale the average too much.
		long double stddevs[sumCount];
		std::vector<long double> stddevForMedian[polarizationCount];
		for(size_t p=0; p!=polarizationCount; ++p)
			stddevForMedian[p].resize(channelsPerSubband);
		size_t sumIndex = 0;
		for(size_t channelIndex=0; channelIndex!=channelsPerSubband; ++channelIndex)
		{
			for(size_t p=0;p!=polarizationCount;++p)
			{
				long double
					rval = stddev(sums[sumIndex].real(), sumsSq[sumIndex].real(), counts[sumIndex]),
					ival = stddev(sums[sumIndex].imag(), sumsSq[sumIndex].imag(), counts[sumIndex]);
				long double avgval = 0.5 * (rval + ival);
				stddevs[sumIndex] = avgval;
				stddevForMedian[p][channelIndex] = avgval;
				++sumIndex;
			}
		}
		
		long double normFactors[polarizationCount];
		if(doNormalize)
		{
			size_t midElement = channelsPerSubband/2;
			for(size_t p=0;p!=polarizationCount;++p)
			{
				std::nth_element(stddevForMedian[p].begin(), stddevForMedian[p].begin()+midElement, stddevForMedian[p].end());
				normFactors[p] = 1.0 / stddevForMedian[p][midElement];
			}
		} else {
			for(size_t p=0;p!=polarizationCount;++p) normFactors[p] = 1.0;
		}
		sumIndex = 0;
		std::cout.precision(12);
		for(size_t channelIndex=0; channelIndex!=channelsPerSubband; ++channelIndex)
		{
			std::cout << channelIndex;
			for(size_t polarizationIndex=0;polarizationIndex!=polarizationCount;++polarizationIndex)
			{
				std::cout
					 << '\t' << (stddevs[sumIndex] * normFactors[polarizationIndex]);
				++sumIndex;
			}
			std::cout << '\n';
		}
	}
	
	return 0;
}
