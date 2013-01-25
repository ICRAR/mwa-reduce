#include <iostream>

#include <string.h>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <stdexcept>

#include <cmath>
#include <fstream>
#include "banddata.h"
#include "predicter.h"
#include "model.h"

using namespace casa;

typedef std::complex<long double> lcomplex_t;
typedef std::complex<float> complex_t;

struct BaselineData {
		BaselineData(size_t timestepCount, size_t channelCount) : _channelCount(channelCount)
		{
			_dataX = new complex_t[timestepCount*channelCount];
			_dataY = new complex_t[timestepCount*channelCount];
		}
		~BaselineData()
		{
			delete[] _dataX;
			delete[] _dataY;
		}
		
		complex_t &DataX(size_t timeIndex, size_t freqIndex) { return _dataX[timeIndex*_channelCount + freqIndex]; }
		complex_t &DataY(size_t timeIndex, size_t freqIndex) { return _dataX[timeIndex*_channelCount + freqIndex]; }
		const complex_t &DataX(size_t timeIndex, size_t freqIndex) const { return _dataX[timeIndex*_channelCount + freqIndex]; }
		const complex_t &DataY(size_t timeIndex, size_t freqIndex) const { return _dataX[timeIndex*_channelCount + freqIndex]; }
		
	private:
		complex_t *_dataX, *_dataY;
		size_t _channelCount;
};

struct BaselineWeights {
		BaselineWeights(size_t timestepCount, size_t channelCount) : _channelCount(channelCount)
		{
			_weightX = new float[timestepCount*channelCount];
			_weightY = new float[timestepCount*channelCount];
			for(size_t i=0; i!=timestepCount*channelCount; ++i)
			{
				_weightX[i] = 0.0;
				_weightY[i] = 0.0;
			}
		}
		~BaselineWeights()
		{
			delete[] _weightX;
			delete[] _weightY;
		}
				
		float &WeightX(size_t timeIndex, size_t freqIndex) { return _weightX[timeIndex*_channelCount + freqIndex]; }
		float &WeightY(size_t timeIndex, size_t freqIndex) { return _weightY[timeIndex*_channelCount + freqIndex]; }
		const float &WeightX(size_t timeIndex, size_t freqIndex) const { return _weightX[timeIndex*_channelCount + freqIndex]; }
		const float &WeightY(size_t timeIndex, size_t freqIndex) const { return _weightY[timeIndex*_channelCount + freqIndex]; }
	private:
		size_t _channelCount;
		float *_weightX, *_weightY;
};

struct Data {
	private:
		BaselineData **_predicted, **_measured;
		BaselineWeights **_weights;
		size_t _antennaCount, _timestepCount, _channelCount;
	
	public:
		Data(size_t antennaCount, size_t timestepCount, size_t channelCount) :
		_antennaCount(antennaCount), _timestepCount(timestepCount), _channelCount(channelCount)
		{
			_predicted = new BaselineData*[antennaCount*antennaCount];
			_measured = new BaselineData*[antennaCount*antennaCount];
			_weights = new BaselineWeights*[antennaCount*antennaCount];
			for(size_t a1=0; a1!=antennaCount; ++a1)
			{
				for(size_t a2=a1+1; a2!=antennaCount; ++a2)
				{
					_predicted[a1*_antennaCount + a2] = new BaselineData(timestepCount, channelCount);
					_measured[a1*_antennaCount + a2] = new BaselineData(timestepCount, channelCount);
					_weights[a1*_antennaCount + a2] = new BaselineWeights(timestepCount, channelCount);
				}
			}
		}
		~Data()
		{
			for(size_t a1=0; a1!=_antennaCount; ++a1)
			{
				for(size_t a2=a1+1; a2!=_antennaCount; ++a2)
				{
					delete &Predicted(a1, a2);
					delete &Measured(a1, a2);
					delete &Weights(a1, a2);
				}
			}
			delete[] _predicted;
			delete[] _measured;
		}
		size_t AntennaCount() const { return _antennaCount; }
		size_t TimestepCount() const { return _timestepCount; }
		size_t ChannelCount() const { return _channelCount; }
		BaselineData &Predicted(size_t a1, size_t a2) { return *_predicted[a1*_antennaCount + a2]; }
		BaselineData &Measured(size_t a1, size_t a2) { return *_measured[a1*_antennaCount + a2]; }
		BaselineWeights &Weights(size_t a1, size_t a2) { return *_weights[a1*_antennaCount + a2]; }
		const BaselineData &Predicted(size_t a1, size_t a2) const { return *_predicted[a1*_antennaCount + a2]; }
		const BaselineData &Measured(size_t a1, size_t a2) const { return *_measured[a1*_antennaCount + a2]; }
		const BaselineWeights &Weights(size_t a1, size_t a2) const { return *_weights[a1*_antennaCount + a2]; }
};

/**
 * Calculate difference between phases and make sure it's
 * within -pi ... pi .
 */
long double PhaseOffset(long double base, long double other)
{
	long double offset = base - other;
	while(offset > M_PIl) offset -= 2.0*M_PIl;
	while(offset <= -M_PIl) offset += 2.0*M_PIl;
	return offset;
}

/**
 * Calculate difference between phases and make sure it's
 * within -pi ... pi .
 */
long double PhaseDistUp(long double base, long double other)
{
	long double offset = base - other;
	while(offset >= 2.0*M_PIl) offset -= 2.0*M_PIl;
	while(offset < 0.0) offset += 2.0*M_PIl;
	return offset;
}

std::string PhaseToString(long double phase)
{
	std::stringstream s;
	double val = ((phase/M_PIl)*180.0L);
	if(val < 0.0)
	{
		s << '-';
		val = -val;
	}
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

lcomplex_t rotate(complex_t data, long double rotSin, long double rotCos)
{
	return lcomplex_t(
		data.real() * rotCos - data.imag() * rotSin,
		data.real() * rotSin + data.imag() * rotCos);
}

void CalculateTimeIntegratedPhaseDifferences(long double **phasesPerBaseline, long double **visWeights, const Data &dataSet, const long double *const *phasesPerAntenna)
{
	const size_t
		antennaCount = dataSet.AntennaCount(),
		channelCount = dataSet.ChannelCount(),
		timestepCount = dataSet.TimestepCount();
	for(size_t a1=0;a1!=antennaCount;++a1)
	{
		for(size_t a2=a1+1;a2!=antennaCount;++a2)
		{
			const long double *antenna1Phases = phasesPerAntenna[a1];
			const long double *antenna2Phases = phasesPerAntenna[a2];
			const BaselineData &mData = dataSet.Measured(a1, a2);
			const BaselineData &pData = dataSet.Predicted(a1, a2);
			const BaselineWeights &wData = dataSet.Weights(a1, a2);
			
			size_t index = SumIndexA2highest(a1, a2, antennaCount);
			
			long double *baselinePhases = phasesPerBaseline[index];
			long double *baselineWeights = visWeights[index];
			
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				long double diffSumXa = 0.0, diffSumYa = 0.0;
				long double diffSumXb = 0.0, diffSumYb = 0.0;
				
				long double curPhaseSolutionX = *antenna1Phases - *antenna2Phases;
				long double curPhaseSolutionY = *(antenna1Phases+1) - *(antenna2Phases+1);
				long double &baselPhaseX = *baselinePhases, &baselPhaseY = *(baselinePhases+1);
				long double wX = 0.0, wY = 0.0;
				
				for(size_t t=0; t!=timestepCount; ++t)
				{
					complex_t mX = mData.DataX(t, ch), mY = mData.DataY(t, ch);
					complex_t pX = pData.DataX(t, ch), pY = pData.DataY(t, ch);
					long double
						phaseX = PhaseOffset(atan2l(mX.imag(), mX.real()), curPhaseSolutionX),
						phaseY = PhaseOffset(atan2l(mY.imag(), mY.real()), curPhaseSolutionY),
						predPhaseX = atan2l(pX.imag(), pX.real()),
						predPhaseY = atan2l(pY.imag(), pY.real());
					long double weightX = wData.WeightX(t, ch);
					long double weightY = wData.WeightY(t, ch);
					long double diffX = PhaseOffset(predPhaseX, phaseX);
					long double diffY = PhaseOffset(predPhaseY, phaseY);
					diffSumXa += diffX * weightX;
					diffSumYa += diffY * weightY;
					diffSumXb += (2.0*M_PIl-diffX) * weightX;
					diffSumYb += (2.0*M_PIl-diffY) * weightY;
					wX += weightX;
					wY += weightY;
				}
				long double diffSumX, diffSumY;
				
				//if(diffSumXa < diffSumXb)
					diffSumX = diffSumXa;
				//else
				//	diffSumX = -diffSumXb;
				
				//if(diffSumYa < diffSumYb)
					diffSumY = diffSumYa;
				//else
				//	diffSumY = -diffSumYb;
				
				if(wX != 0.0)
					baselPhaseX = diffSumX / wX;
				else
					baselPhaseX = 0.0;
				if(wY != 0.0)
					baselPhaseY = diffSumY / wY;
				else
					baselPhaseY = 0.0;
				++baselinePhases;
				++baselinePhases;
				++baselineWeights;
				*baselineWeights = wX;
				++baselineWeights;
				*baselineWeights = wY;
				antenna1Phases += 2;
				antenna2Phases += 2;
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

long double AntennaPhaseOffset(size_t antenna, size_t antennaCount, size_t eIndex, const long double* const* phasesPerBaseline, const long double* const* phaseWeights)
{
	long double phaseOffset = 0, weightSum = 0.0;
	// a2 < a1 ==> a1 is not conjugated.
	for(size_t a2 = 0; a2 != antenna; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		const long double *baselineValues = phasesPerBaseline[aIndex];
		const long double weight = phaseWeights[aIndex][eIndex];
		const long double posPhase = baselineValues[eIndex];
		phaseOffset -= posPhase * weight;
		weightSum += weight;
	}
	// a2 > a1 ==> a1 is conjugated.
	for(size_t a2 = antenna+1; a2!=antennaCount; ++a2)
	{
		size_t aIndex = SumIndex(antenna, a2, antennaCount);
		const long double *baselineValues = phasesPerBaseline[aIndex];
		const long double weight = phaseWeights[aIndex][eIndex];
		const long double posPhase = PhaseDistUp(-baselineValues[eIndex], 0.0);
		phaseOffset += posPhase * weight;
		weightSum += weight;
	}
	if(antenna==6 && eIndex==2)
	{
		//std::cout << PhaseToString(phaseOffsetA / weightSum) << "," << PhaseToString(-phaseOffsetB / weightSum) << '\n';
	}
	if(weightSum != 0.0)
	{
		return phaseOffset / weightSum;
	}
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

std::pair<long double, long double> ChannelStdError(size_t ch, size_t antennaCount, size_t polarizationCount, const long double* const* phasesPerBaseline, const long double* const* phaseWeights)
{
	long double stdError = 0.0, weightSum = 0.0;
	size_t eIndex = polarizationCount * ch;
	for(size_t pol=0; pol!=polarizationCount; ++pol)
	{
		for(size_t a1=0;a1!=antennaCount;++a1)
		{
			long double phaseOffset = AntennaPhaseOffset(a1, antennaCount, eIndex, phasesPerBaseline, phaseWeights);
			
			long double weight = AntennaWeight(a1, antennaCount, eIndex, phaseWeights);
			stdError += phaseOffset * phaseOffset * weight;
			weightSum += weight;
		}
	}
	if(weightSum == 0.0)
		return std::pair<long double, long double>(0.0, 0.0);
	else
		return std::pair<long double, long double>(sqrtl(stdError / weightSum), weightSum);
}

void FitPhases(size_t channelCount, size_t polarizationCount, long double *phaseOffsets, long double *outputs, size_t outputChannelCount)
{
	size_t eTotal = channelCount * polarizationCount;
	for(size_t p=0;p!=polarizationCount;++p)
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

void save(const char *outName, size_t antennaCount, size_t avgChannelCount, size_t inpChannelCount, size_t polarizationCount, const long double* const* phases, bool fitSlope)
{
	std::ofstream outFile(outName);
	outFile.precision(10);
	outFile << antennaCount << '\t';
	if(avgChannelCount == 1)
	{
		outFile << avgChannelCount << '\t' << polarizationCount << '\n';
		for(size_t a = 0; a!=antennaCount; ++a)
		{
			outFile << a;
			for(size_t p = 0; p!=polarizationCount; ++p)
			{
				if(p == 0 || p == polarizationCount-1)
				{
					size_t pIndex = (p==0) ? 0 : 1;
					const long double *antennaPhases = phases[a];
					outFile << '\t' << antennaPhases[pIndex];
				} else {
					outFile << "\t0.0";
				}
			}
			outFile << '\n';
		}
	} else {
		size_t eIndex = 0;
		size_t outpChannelCount = fitSlope ? inpChannelCount : avgChannelCount;
		outFile << outpChannelCount << '\t' << polarizationCount << '\n';
		for(size_t ch = 0; ch!=outpChannelCount; ++ch)
		{
			outFile << ch;
			for(size_t p = 0; p!=polarizationCount; ++p)
			{
				if(p == 0 || p == polarizationCount-1)
				{
					for(size_t a = 0; a!=antennaCount; ++a)
					{
						const long double *antennaPhases = phases[a];
						outFile << '\t' << antennaPhases[eIndex];
					}
					++eIndex;
				} else {
					for(size_t a = 0; a!=antennaCount; ++a)
						outFile << "\t0.0";
				}
			}
			outFile << '\n';
		}
	}
}

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: calcphaseoffsets [-a <nr>] [-f] <model> <measurementset.ms> <output.txt>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n"
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
		if(argc <= argi + 2) throw std::runtime_error("Incorrect parameters");
		const char *modelName = argv[argi];
		const char *msName = argv[argi+1];
		const char *outName = argv[argi+2];
		MeasurementSet ms(msName);
		
		std::cout << "Reading model... " << std::flush;
		Model model(modelName);
		std::cout << "DONE\n";
		
		std::cout << "Reading measurement set... " << std::flush;
		
		/**
		 * Read some meta data from the measurement set
		 */
		MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		if(avgChannelCount == 0) avgChannelCount = channelCount;
		
		MSField fieldTable = ms.field();
		ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ROScalarColumn<double> timeColumn(ms, ms.columnName(MSMainEnums::TIME));
		ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ROArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		ROArrayColumn<double> uvwColumn(ms, ms.columnName(MSMainEnums::UVW));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		std::cout << "DONE\nCounting timesteps... " << std::flush;
		double time = -1.0;
		size_t timestepCount = 0;
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			if(timeColumn(rowIndex) != time)
			{
				++timestepCount;
				time = timeColumn(rowIndex);
			}
		}
		
		/**
		 * Initialize memory
		 */
		const size_t matrixSize = antennaCount * antennaCount;
		size_t sumsPerBaseline = avgChannelCount * 2;
		size_t reservedSpace = channelCount * 2; // we reserve more for interpolation
		long double *phaseValues[matrixSize], *phaseWeights[matrixSize];
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			phaseValues[e] = new long double[reservedSpace];
			phaseWeights[e] = new long double[reservedSpace];
			for(size_t i=0; i!=reservedSpace; ++i)
			{
				phaseValues[e][i] = 0.0;
				phaseWeights[e][i] = 0.0;
			}
		}
		
		long double *phasesPerAntenna[antennaCount];
		for(size_t a=0; a!=antennaCount; ++a) 
		{
			phasesPerAntenna[a] = new long double[reservedSpace];
			for(size_t i=0; i!=reservedSpace; ++i) phasesPerAntenna[a][i] = 0.0;
		}
		Data dataSet(antennaCount, timestepCount, avgChannelCount);
		
		std::cout << "DONE (" << timestepCount << ")\nReading data... " << std::flush;
		
		/**
		 * Read the data. Average channels if requested.
		 */
		{
			Array<complex_t> data(dataShape);
			Array<bool> flags(dataShape);
			size_t timeIndex = 0;
			time = timeColumn(0);
			const double sampleWeight = 1.0 / channelCount;
			for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
			{
				if(timeColumn(rowIndex) != time)
				{
					++timeIndex;
					time = timeColumn(rowIndex);
				}
				// Cross correlation?
				size_t antenna1 = ant1Column.get(rowIndex), antenna2 = ant2Column.get(rowIndex);
				if(antenna1 != antenna2)
				{
					dataColumn.get(rowIndex, data);
					flagColumn.get(rowIndex, flags);
					Array<complex_t>::const_iterator dataPtr = data.begin();
					Array<bool>::const_iterator flagPtr = flags.begin();
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						size_t chIndex = (ch*avgChannelCount/channelCount);
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							// Only include non-flagged data
							if(!(*flagPtr) && std::isfinite(dataPtr->real()) && std::isfinite(dataPtr->imag())) {
								if(p == 0)
								{
									dataSet.Measured(antenna1, antenna2).DataX(timeIndex, chIndex) += *dataPtr;
									dataSet.Weights(antenna1, antenna2).WeightX(timeIndex, chIndex) += sampleWeight;
								}
								else if(p == polarizationCount-1)
								{
									dataSet.Measured(antenna1, antenna2).DataY(timeIndex, chIndex) += *dataPtr;
									dataSet.Weights(antenna1, antenna2).WeightY(timeIndex, chIndex) += sampleWeight;
								}
							}
							++dataPtr;
							++flagPtr;
						}
					}
				}
			}
		}
		
		std::cout << "DONE!\nPredicting model... " << std::flush;
		Predicter predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), avgChannelCount);
		predicter.Initialize(model);
		size_t timeIndex = 0;
		time = timeColumn(0);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			size_t antenna1 = ant1Column.get(rowIndex), antenna2 = ant2Column.get(rowIndex);
			if(timeColumn(rowIndex) != time)
			{
				++timeIndex;
				time = timeColumn(rowIndex);
			}
			if(antenna1 != antenna2)
			{
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				BaselineData &predictedData = dataSet.Predicted(antenna1, antenna2);
				for(size_t ch = 0; ch!=avgChannelCount; ++ch)
				{
					double lambda = bandData.ChannelWavelength(ch*channelCount/avgChannelCount);
					lcomplex_t p = predicter.Predict(model, u/lambda, v/lambda, w/lambda, ch);
					predictedData.DataX(timeIndex, ch) = p;
					predictedData.DataY(timeIndex, ch) = p;
				}
			}
		}
		
		/**
		 *  Calculate average phase per baseline per channel (differential phase of average vis)
		 */
		
		std::cout << "DONE!\nInitial standard error of phase offsets: " << std::flush;
		
		size_t inpChannelCount = channelCount;

		CalculateTimeIntegratedPhaseDifferences(phaseValues, phaseWeights, dataSet, phasesPerAntenna);

		long double stdError = 0.0, stdErrorWeight = 0.0;
		for(size_t ch=0; ch!=avgChannelCount; ++ch)
		{
				std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, 2, phaseValues, phaseWeights);
				stdError += result.first * result.second;
				stdErrorWeight += result.second;
		}
		stdError /= stdErrorWeight;
		std::cout << PhaseToString(stdError) << ".\n";
		
		size_t iteration=0, gotWorseCount=0;
		long double prevAbsError = 10.0, minError = 10.0, stepSize = 0.1;
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
				for(size_t ch = 0; ch!=avgChannelCount; ++ch)
				{
					for(size_t p = 0; p!=2; ++p)
					{
						phaseOffsets[eIndex] = AntennaPhaseOffset(a1, antennaCount, eIndex, phaseValues, phaseWeights);
						antennaWeights[eIndex] = AntennaWeight(a1, antennaCount, eIndex, phaseWeights);
						
						if(eIndex/2+1 == 3 && a1==6)
						{
							std::cout << PhaseToString(phaseOffsets[eIndex]) << '\t' << PhaseToString(phasesPerAntenna[a1][eIndex]) << '\t' << "W=" << antennaWeights[eIndex] << '\t';
						}
						++eIndex;
					}
				}
				
				// Apply the average phase offset to the 'phaseValues'
				eIndex = 0;
				for(size_t ch = 0; ch!=avgChannelCount; ++ch)
				{
					for(size_t p = 0; p!=2; ++p)
					{
						long double step = phaseOffsets[eIndex] * stepSize * antennaWeights[eIndex];
						phasesPerAntenna[a1][eIndex] += step;
						++eIndex;
					}
				}
			}
			
			CalculateTimeIntegratedPhaseDifferences(phaseValues, phaseWeights, dataSet, phasesPerAntenna);
			
			for(size_t ch=0; ch!=avgChannelCount; ++ch)
			{
				std::pair<long double, long double> result = ChannelStdError(ch, antennaCount, 2, phaseValues, phaseWeights);
				stdError += result.first * result.second;
				stdErrorWeight += result.second;
			}
			stdError /= stdErrorWeight;
			if(stdError < minError) minError = stdError;
			std::cout << "Standard error: " << PhaseToString(stdError) << '\n';
			
			bool gotBetter = stdError < prevAbsError;
			if(!gotBetter) { ++gotWorseCount; stepSize*=0.9; }
			++iteration;

			stringstream solNameStr;
			solNameStr << "temp/sol-";
			if(iteration < 10) solNameStr << '0';
			if(iteration < 100) solNameStr << '0';
			if(iteration < 1000) solNameStr << '0';
			solNameStr << iteration << ".txt";
			save(solNameStr.str().c_str(), antennaCount, avgChannelCount, inpChannelCount, polarizationCount, phasesPerAntenna, false);
			
			stringstream offsetNameStr;
			offsetNameStr << "temp/offset-";
			if(iteration < 10) offsetNameStr << '0';
			if(iteration < 100) offsetNameStr << '0';
			if(iteration < 1000) offsetNameStr << '0';
			offsetNameStr << iteration << ".txt";
			//save(tempNameStr.str().c_str(), antennaCount, avgChannelCount, inpChannelCount, polarizationCount, phaseOffsets, false);
		} while(iteration < 250); // stdError > 0.00001  && 
		
		if(fitSlope) {
			for(size_t a = 0; a!=antennaCount; ++a)
			{
				long double *antennaPhases = phasesPerAntenna[a];
				FitPhases(avgChannelCount, 2, antennaPhases, antennaPhases, inpChannelCount);
			}
		}
		
		save(outName, antennaCount, avgChannelCount, inpChannelCount, polarizationCount, phasesPerAntenna, fitSlope);
				
		for(size_t e=0; e!=matrixSize; ++e) 
		{
			delete[] phaseValues[e];
			delete[] phaseWeights[e];
		}
		for(size_t a=0; a!=antennaCount;++a)
			delete[] phasesPerAntenna[a];
	}
	
	return 0;
}
