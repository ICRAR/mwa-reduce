#include "msprovider.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "../msselection.h"

void MSProvider::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(polOut == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			} else {
				dest[ch] = 0;
			}
			weightPtr += polCount-1;
			inPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] += *inPtr * (*weightPtr);
			}
			++weightPtr;
			++inPtr;
			++flagPtr;
		}
	} /*else if(Polarization() == Polarization::XY || Polarization() == Polarization::YX)
	{
		// Step to XY:
		++weightPtr;
		++inPtr;
		++flagPtr;
		const bool flipSign = Polarization() == Polarization::YX;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			// Step to YX:
			++weightPtr;
			++inPtr;
			++flagPtr;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				//dest[ch] += *inPtr * (*weightPtr) * rowWeight;
				//_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 3;
			inPtr += 3;
			flagPtr += 3;
			//if(flipSign)
			//	dest[ch].imag(dest[ch].imag() * -1.0);
		}
	}*/ else {
		int polIndex = Polarization::TypeToIndex(polOut, polCount);
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr);
			}
			else {
				dest[ch] = 0;
			}
			weightPtr += polCount;
			inPtr += polCount;
			flagPtr += polCount;
		}
	}
}

template<typename NumType>
void MSProvider::copyWeights(NumType* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(polOut == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polCount-1;
			weightPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] += *weightPtr;
			++inPtr;
			++weightPtr;
			++flagPtr;
		}
	} /*else if(Polarization() == Polarization::XY || Polarization() == Polarization::YX)
	{
		// Step to XY:
		inPtr++;
		weightPtr++;
		flagPtr++;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			// Step to YX:
			inPtr++;
			weightPtr++;
			flagPtr++;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				//dest[ch] += (*weightPtr) * rowWeight;
				//_totalWeight += (*weightPtr) * rowWeight;
			}
			inPtr += 3;
			weightPtr += 3;
			flagPtr += 3;
		}
	} */ else {
		int polIndex = Polarization::TypeToIndex(polOut, polCount);
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polCount;
			weightPtr += polCount;
			flagPtr += polCount;
		}
	}
}

template
void MSProvider::copyWeights<float>(float* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

template
void MSProvider::copyWeights<std::complex<float>>(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

void MSProvider::reverseCopyData(casa::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, size_t polCount, const std::complex<float>* source, PolarizationEnum polSource)
{
	const size_t selectedChannelCount = endChannel - startChannel;
	casa::Array<std::complex<float>>::contiter dataIter = dest.cbegin() + startChannel * polCount;
	
	if(polSource == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*dataIter = source[ch];
				*(dataIter + (polCount-1)) = source[ch];
			}
			dataIter += polCount;
		}
	} else {
		int polIndex = Polarization::TypeToIndex(polSource, polCount);
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*(dataIter+polIndex) = source[ch];
			}
			dataIter += polCount;
		}
	}
}

void MSProvider::getRowRange(casa::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow)
{
	startRow = 0;
	endRow = ms.nrow();
	if(selection.HasInterval())
	{
		std::cout << "Determining first and last row index... " << std::flush;
		casa::ROScalarColumn<double> timeColumn(ms, casa::MS::columnName(casa::MSMainEnums::TIME));
		double time = timeColumn(0);
		size_t timestepIndex = 0;
		for(size_t row = 0; row!=ms.nrow(); ++row)
		{
			if(time != timeColumn(row))
			{
				++timestepIndex;
				if(timestepIndex == selection.IntervalStart())
					startRow = row;
				if(timestepIndex == selection.IntervalEnd())
				{
					endRow = row;
					break;
				}
				time = timeColumn(row);
			}
		}
		std::cout << "DONE (" << startRow << '-' << endRow << ")\n";
	}
}

