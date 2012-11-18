#include <iostream>

#include <string.h>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <stdexcept>

#include <cmath>
#include <fstream>

using namespace casa;

typedef std::complex<long double> lcomplex_t;
typedef std::complex<float> complex_t;

/**
 * Calculate difference between phases and make sure it's
 * within -pi ... pi .
 */
long double PhaseOffset(long double base, long double other)
{
	long double offset = base - other;
	if(offset > M_PIl) offset -= 2.0*M_PIl;
	if(offset <= -M_PIl) offset += 2.0*M_PIl;
	return offset;
}

std::string PhaseToString(long double phase)
{
	std::stringstream s;
	double val = ((phase/M_PIl)*180.0L);
	const char *unit;
	if(val >= 1.0)
	{
		unit = "deg";
	} else if(val*60.0>=1.0)
	{
		val *= 60.0;
		unit = "amin";
	} else if(val*60*60>=1.0)
	{
		val *= 60.0*60.0;
		unit = "asec";
	} else {
		val *= 60.0*60.0*1000.0;
		unit = "masec";
	}
	s.precision(3);
	s << val << ' ' << unit;
	return s.str();
}

size_t SumIndex(size_t a1, size_t a2, size_t antennaCount)
{
	return (a1>a2) ? (a2*antennaCount + a1) : (a1*antennaCount + a2);
}
size_t SumIndexA2highest(size_t a1, size_t a2, size_t antennaCount)
{
	return a1*antennaCount + a2;
}
size_t SumIndexA1highest(size_t a1, size_t a2, size_t antennaCount)
{
	return a2*antennaCount + a1;
}

void SumsToPhases(long double **phases, const lcomplex_t* const* sums, size_t antennaCount, size_t polarizationCount, size_t channelCount)
{
	for(size_t a1=0;a1!=antennaCount;++a1)
	{
		for(size_t a2=a1+1;a2!=antennaCount;++a2)
		{
			size_t index = SumIndexA2highest(a1, a2, antennaCount);
			long double *antPhases = phases[index];
			
			const lcomplex_t *antSums = sums[index];
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					*antPhases = atan2l(antSums->imag(), antSums->real());
					++antPhases;
					++antSums;
				}
			}
		}
	}
}

void UnwrapPhase(long double *phases, size_t count, size_t step)
{
	const long double twopi = 2.0*M_PIl;
	long double prevPhase = 0.0;
	for(size_t i=0;i!=count;++i)
	{
		*phases = fmodl(*phases-prevPhase, twopi)+prevPhase;
		if(prevPhase - *phases > M_PIl)
			*phases += twopi;
		else if(*phases - prevPhase > M_PIl)
			*phases -= twopi;
		if(std::isfinite(*phases))
			prevPhase = *phases;
		phases += step;
	}
}

long double AntennaPhaseOffset(size_t antenna, size_t antennaCount, size_t eIndex, const long double* const* phaseValues, const long double* const* phaseWeights, long double refPhase = 0.0)
{
	long double phaseOffset = 0, weightSum = 0.0;
	// a2 < a1 ==> a1 is not conjugated.
	for(size_t a2 = 0; a2 != antenna; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		const long double *baselineValues = phaseValues[aIndex];
		const long double weight = phaseWeights[aIndex][eIndex];
		phaseOffset += PhaseOffset(refPhase, baselineValues[eIndex]) * weight;
		weightSum += weight;
	}
	// a2 > a1 ==> a1 is conjugated.
	for(size_t a2 = antenna+1; a2!=antennaCount; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		const long double *baselineValues = phaseValues[aIndex];
		const long double weight = phaseWeights[aIndex][eIndex];
		phaseOffset -= PhaseOffset(refPhase, baselineValues[eIndex]) * weight;
		weightSum += weight;
	}
	if(weightSum != 0.0)
		return phaseOffset / weightSum;
	else
		return 0.0;
}

long double AntennaWeight(size_t antenna, size_t antennaCount, size_t eIndex, const long double* const* visWeights)
{
	long double weight = 0;
	for(size_t a2 = 0; a2 != antenna; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		weight += visWeights[aIndex][eIndex];
	}
	for(size_t a2 = antenna+1; a2 != antennaCount; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		weight += visWeights[aIndex][eIndex];
	}
	weight /= (long double) (antennaCount-1);
	return weight;
}

std::pair<long double, long double> ChannelStdError(size_t ch, size_t antennaCount, size_t polarizationCount, const long double* const* phaseValues, const long double* const* phaseWeights)
{
	long double stdError = 0.0, weightSum = 0.0;
	size_t eIndex = polarizationCount * ch;
	for(size_t pol=0; pol!=polarizationCount; ++pol)
	{
		if(polarizationCount!=4 || pol==0 || pol==3) {
			for(size_t a1=0;a1!=antennaCount;++a1)
			{
				long double phaseOffset = AntennaPhaseOffset(a1, antennaCount, eIndex, phaseValues, phaseWeights, 0.0);
				
				if(polarizationCount!=4 || pol==0 || pol==3)
				{
					long double weight = AntennaWeight(a1, antennaCount, eIndex, phaseWeights);
					stdError += phaseOffset * phaseOffset * weight;
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

long double TotalPhaseOffset(const lcomplex_t* const* sums, size_t antennaCount, size_t eIndex)
{
	/* Function should take into account conjugate! But function is not used, so
	 * leave it for now. */
	lcomplex_t antennaSum;
	for(size_t a1 = 0; a1!=antennaCount; ++a1)
	{
		for(size_t a2 = 0; a2!=antennaCount; ++a2)
		{
			size_t aIndex = SumIndex(a1, a2, antennaCount);
			antennaSum += sums[aIndex][eIndex];
		}
	}
	// Calibrate towards refPhase. If the source is not in phase centre, we could
	// use this:
	// long double refPhase = atan2l(antennaSum.imag(), antennaSum.real());
	// (but is much less accurate and stable)
	return atan2(antennaSum.imag(), antennaSum.real());
}

void FitPhases(size_t channelCount, size_t polarizationCount, long double *phaseOffsets, long double *outputs, size_t outputChannelCount)
{
	size_t eTotal = channelCount * polarizationCount;
	for(size_t p=0;p!=polarizationCount;++p)
	{
		if(p!=polarizationCount || (p==0 || p==3))
		{
			UnwrapPhase(phaseOffsets+p, channelCount, polarizationCount);
			
			long double avgY = 0.0, avgX = (long double) (channelCount-1.0) / 2.0;
			for(size_t eIndex=p;eIndex<eTotal;eIndex+=polarizationCount)
				avgY += phaseOffsets[eIndex];
			avgY /= (long double) channelCount;
			
			long double ssxx = 0.0, ssxy = 0.0;
			size_t x = 0;
			for(size_t eIndex=p;eIndex<eTotal;eIndex+=polarizationCount)
			{
				long double
					diffx = (long double) x - avgX,
					diffy = phaseOffsets[eIndex] - avgY;
				ssxx += diffx*diffx;
				ssxy += diffx*diffy;
				++x;
			}
			long double beta = (ssxx!=0.0) ? ssxy / ssxx : 0.0;
			long double alpha = avgY - beta * avgX;
			x = 0;
			long double xFactor = beta * (long double) channelCount / outputChannelCount;
			long double xSum = ((long double) outputChannelCount / channelCount)/2.0L - 0.5L;
			xSum = xSum*xFactor + alpha;
			for(size_t eIndex=p;eIndex<outputChannelCount*polarizationCount;eIndex+=polarizationCount)
			{
				outputs[eIndex] = (long double) x * xFactor + xSum;
				++x;
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
			size_t index = SumIndexA2highest(a1, a2, antennaCount);
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
			size_t index = SumIndexA2highest(a1, a2, antennaCount);
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

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		std::cout
			<< "Usage: calcphaseoffsets [-a <nr>] <measurementset.ms> <output.txt>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. To do this, it assumes the array\n"
			<< "is tracking a bright source in the phase centre. It produces approximate least-squares solutions.\n"
			<< "Option -a will average over frequency before fitting, nr should specify the amount\n"
			<< "of desired channels.\n";
	} else {
		bool fitSlope = false;
		int argi = 1;
		size_t avgChannelCount = 0;
		while(argi<argc && argv[argi][0]=='-')
		{
			if(strcmp(argv[argi], "-a")==0) {
				++argi;
				avgChannelCount = atoi(argv[argi]);
				++argi;
			}
			else if(strcmp(argv[argi], "-f")==0) {
				fitSlope = true;
				++argi;
			}
		}
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		MeasurementSet ms(msName);
		
		std::cout << "Reading measurement set... " << std::flush;
		
		/**
		 * Read some meta data from the measurement set
		 */
		MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		MSSpectralWindow spwTable = ms.spectralWindow();
		size_t spwCount = spwTable.nrow();
		if(spwCount != 1) throw std::runtime_error("Set should have exactly one spectral window");
		
		ROScalarColumn<int> numChanCol(spwTable, MSSpectralWindow::columnName(MSSpectralWindowEnums::NUM_CHAN));
		size_t channelCount = numChanCol.get(0);
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		if(avgChannelCount == 0) avgChannelCount = channelCount;
		
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
		lcomplex_t *sums[matrixSize];
		long double *phaseValues[matrixSize], *phaseWeights[matrixSize];
		size_t *counts[matrixSize];
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			sums[e] = new lcomplex_t[sumsPerBaseline];
			phaseValues[e] = new long double[sumsPerBaseline];
			phaseWeights[e] = new long double[sumsPerBaseline];
			counts[e] = new size_t[sumsPerBaseline];
			memset(&counts[e][0], 0, sizeof(size_t) * sumsPerBaseline);
		}
		
		long double *phases[antennaCount];
		for(size_t a=0; a!=antennaCount; ++a) 
		{
			phases[a] = new long double[sumsPerBaseline];
			for(size_t i=0; i!=sumsPerBaseline; ++i) phases[a][i] = 0.0;
		}
		
		/**
		 * Read the average phase per baseline per channel (phase of average vis)
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
					long double conjSign = 1.0;
					if(antenna1 > antenna2) { std::swap(antenna1, antenna2); conjSign = -1.0; }
					dataColumn.get(rowIndex, data);
					flagColumn.get(rowIndex, flags);
					Array<complex_t>::const_iterator dataPtr = data.begin();
					Array<bool>::const_iterator flagPtr = flags.begin();
					lcomplex_t *baselineSums = sums[antenna2 + antenna1 * antennaCount];
					size_t *baselineCounts = counts[antenna2 + antenna1 * antennaCount];
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						size_t eIndex = (ch*avgChannelCount/channelCount) * polarizationCount;
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							// Only include non-flagged data
							if(!(*flagPtr) && std::isfinite(dataPtr->real()) && std::isfinite(dataPtr->imag())) {
								long double r = dataPtr->real(), i = dataPtr->imag();
								baselineSums[eIndex] += lcomplex_t(r, i * conjSign);
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
		
		std::cout << "DONE!\nInitial standard error of phase offsets: " << std::flush;
		
		size_t maxCount = MaximumCount(counts, antennaCount, polarizationCount, channelCount);
		InitializeWeights(phaseWeights, counts, antennaCount, polarizationCount, avgChannelCount, maxCount);
		
		size_t inpChannelCount = channelCount;
		channelCount = avgChannelCount;
		sumsPerBaseline = polarizationCount * channelCount;

		SumsToPhases(phaseValues, sums, antennaCount, polarizationCount, channelCount);

		long double stdError = 0.0, stdErrorWeight = 0.0;
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
				std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, polarizationCount, phaseValues, phaseWeights);
				stdError += result.first;
				stdErrorWeight += result.second;
		}
		stdError /= stdErrorWeight;
		std::cout << PhaseToString(stdError) << ".\n";
		
		size_t iteration=0, gotWorseCount=0;
		long double prevAbsError = 10.0, minError = 10.0, stepSize = 1.0;
		do
		{
			std::cout << "Iteration " << iteration << ": " << std::flush;
			prevAbsError = stdError;
			stdError = 0.0;
			stdErrorWeight = 0.0;
			
			for(size_t a1 = 0; a1!=antennaCount; ++a1)
			{
				/**
				* Calculate per antenna phase offsets
				*/
				size_t eIndex = 0;
				long double phaseOffsets[sumsPerBaseline], antennaWeights[sumsPerBaseline];
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					for(size_t p = 0; p!=polarizationCount; ++p)
					{
						// Skip XY / YX if they exist
						if(polarizationCount!=4 || p==0 || p==3)
						{
							long double refPhase = 0.0;
							phaseOffsets[eIndex] = AntennaPhaseOffset(a1, antennaCount, eIndex, phaseValues, phaseWeights, refPhase);
							antennaWeights[eIndex] = AntennaWeight(a1, antennaCount, eIndex, phaseWeights);
							
							if(eIndex/polarizationCount+1 == (channelCount/2) && a1==1)
							{
								std::cout << PhaseToString(phaseOffsets[eIndex]) << '\t' << PhaseToString(phases[a1][eIndex]) << '\t' << "W=" << antennaWeights[eIndex] << '\t';
							}
						}
						++eIndex;
					}
				}
				
				// Apply the average phase offset to the 'phaseValues'
				eIndex = 0;
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					for(size_t p = 0; p!=polarizationCount; ++p)
					{
						// Skip XY / YX if they exist
						if(polarizationCount!=4 || p==0 || p==3)
						{
							long double step = phaseOffsets[eIndex] * stepSize * antennaWeights[eIndex];
							phases[a1][eIndex] += step;
							
							/** a1 > a2 : a1 is not conjugated */
							for(size_t a2 = 0; a2!=a1; ++a2)
							{
								size_t aIndex = SumIndex(a1, a2, antennaCount);
								
								long double *baselineValues = phaseValues[aIndex];
								baselineValues[eIndex] += step;
							}
							/** a1 < a2 : a1 is conjugated */
							for(size_t a2 = a1+1; a2!=antennaCount; ++a2)
							{
								size_t aIndex = SumIndex(a1, a2, antennaCount);
								
								long double *baselineValues = phaseValues[aIndex];
								baselineValues[eIndex] -= step;
							}
						}
						++eIndex;
					}
				}
			}
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, polarizationCount, phaseValues, phaseWeights);
				stdError += result.first;
				stdErrorWeight += result.second;
			}
			stdError /= stdErrorWeight;
			if(stdError < minError) minError = stdError;
			std::cout << "Standard error: " << PhaseToString(stdError) << '\n';
			
			bool gotBetter = stdError < prevAbsError;
			if(!gotBetter) { ++gotWorseCount; stepSize*=0.9; }
			++iteration;
		} while((gotWorseCount < 10 || stdError > minError*1.1) && iteration < 75);
		
		if(fitSlope) {
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				long double *antennaPhases = phases[a];
				FitPhases(channelCount, polarizationCount, antennaPhases, antennaPhases, inpChannelCount);
			}
		}
				
		/**
		* Print results
		*/
		std::ofstream outFile(outName);
		outFile.precision(10);
		outFile << antennaCount << '\t';
		if(channelCount == 1)
		{
			outFile << channelCount << '\t' << polarizationCount << '\n';
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				outFile << a;
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					long double *antennaPhases = phases[a];
					outFile << '\t' << antennaPhases[p];
				}
				outFile << '\n';
			}
		} else {
			size_t eIndex = 0;
			size_t outpChannelCount = fitSlope ? inpChannelCount : channelCount;
			outFile << outpChannelCount << '\t' << polarizationCount << '\n';
			for(size_t ch = 0; ch!=inpChannelCount; ++ch)
			{
				outFile << ch;
				for(size_t p = 0; p!=polarizationCount; ++p)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						long double *antennaPhases = phases[a];
						outFile << '\t' << antennaPhases[eIndex];
					}
					++eIndex;
				}
				outFile << '\n';
			}
		}
		
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			delete[] sums[e];
			delete[] phaseValues[e];
			delete[] counts[e];
		}
		for(size_t a=0; a!=antennaCount;++a)
			delete[] phases[a];
	}
	
	return 0;
}
