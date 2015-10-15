#include "msprovider.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrColDesc.h>

#include "../msselection.h"

void MSProvider::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casacore::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casacore::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casacore::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	size_t polIndex;
	if(Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
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
	else {
		switch(polOut)
		{
		case Polarization::StokesI: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
			bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
			if(!hasXX || !hasYY)
			{
				bool hasRR = Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
				bool hasLL = Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				if(!hasRR || !hasLL)
					throw std::runtime_error("Can not form requested polarization (Stokes I) from available polarizations");
			}
			
			for(size_t ch=0; ch!=selectedChannelCount; ++ch)
			{
				weightPtr += polIndexA;
				inPtr += polIndexA;
				flagPtr += polIndexA;
				
				bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
				casacore::Complex valA = *inPtr * (*weightPtr);
				
				weightPtr += polIndexB - polIndexA;
				inPtr += polIndexB - polIndexA;
				flagPtr += polIndexB - polIndexA;
				
				bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
				if(flagA || flagB)
					dest[ch] = 0.0;
				else
					dest[ch] = (*inPtr * (*weightPtr)) + valA;
				
				weightPtr += polCount - polIndexB;
				inPtr += polCount - polIndexB;
				flagPtr += polCount - polIndexB;
			}
		} break;
		case Polarization::StokesQ: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
			bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
			if(hasXX && hasYY)
			{
				// Convert to StokesQ from XX and YY
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else
						dest[ch] = (*inPtr * (*weightPtr)) - valA;
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesQ from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes Q) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						// Q = RL LR
						dest[ch] = valA *  *inPtr * (*weightPtr);
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		case Polarization::StokesU: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
			bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
			if(hasXY && hasYX)
			{
				// Convert to StokesU from XY and YX
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else
						dest[ch] = valA + (*inPtr * (*weightPtr)); // U = (YX + XY)/2
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesU from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes U) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casacore::Complex diff = (valA - *inPtr * (*weightPtr));
						// U = -i (RL - LR)
						dest[ch] = casacore::Complex(diff.imag(), -diff.real());
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		case Polarization::StokesV: {
			size_t polIndexA=0, polIndexB=0;
			bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
			bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
			if(hasXY && hasYX)
			{
				// Convert to StokesV from XX and YY
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casacore::Complex diff = (*inPtr * (*weightPtr) - valA);
						// V = -i(YX - XY)/2
						dest[ch] = casacore::Complex(diff.imag(), -diff.real());
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
			else {
				// Convert to StokesV from RR and LL
				bool hasRL = Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
				bool hasLR = Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				if(!hasRL || !hasLR)
					throw std::runtime_error("Can not form requested polarization (Stokes V) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casacore::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						// U = RR - LL
						dest[ch] = valA - *inPtr * (*weightPtr);
					}
					
					weightPtr += polCount - polIndexB;
					inPtr += polCount - polIndexB;
					flagPtr += polCount - polIndexB;
				}
			}
		} break;
		default:
			throw std::runtime_error("Could not convert ms polarizations to requested polarization");
		}
	}
}

template<typename NumType>
void MSProvider::copyWeights(NumType* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casacore::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casacore::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casacore::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	size_t polIndex;
	if(Polarization::TypeToIndex(polOut, polsIn, polIndex)) {
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
	else {
		size_t polIndexA=0, polIndexB=0;
		switch(polOut) {
			case Polarization::StokesI: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesQ: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesU: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexB);
				}
			}
			break;
			case Polarization::StokesV: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexB);
				}
			}
			break;
			default:
				throw std::runtime_error("Could not convert ms polarizations to requested polarization");
			break;
		}
		
		weightPtr += polIndexA;
		inPtr += polIndexA;
		flagPtr += polIndexA;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] = *weightPtr;
			else
				dest[ch] = 0;
			inPtr += polIndexB-polIndexA;
			weightPtr += polIndexB-polIndexA;
			flagPtr += polIndexB-polIndexA;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				dest[ch] += *weightPtr;
			weightPtr += polCount - polIndexB + polIndexA;
			inPtr += polCount - polIndexB + polIndexA;
			flagPtr += polCount - polIndexB + polIndexA;
		}
	}
}

template
void MSProvider::copyWeights<float>(float* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);

template
void MSProvider::copyWeights<std::complex<float>>(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casacore::Array<std::complex<float>>& data, const casacore::Array<float>& weights, const casacore::Array<bool>& flags, PolarizationEnum polOut);

void MSProvider::reverseCopyData(casacore::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum> &polsDest, const std::complex<float>* source, PolarizationEnum polSource)
{
	size_t polCount = polsDest.size();
	const size_t selectedChannelCount = endChannel - startChannel;
	casacore::Array<std::complex<float>>::contiter dataIter = dest.cbegin() + startChannel * polCount;
	
	size_t polIndex;
	if(Polarization::TypeToIndex(polSource, polsDest, polIndex)) {
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(std::isfinite(source[ch].real()))
			{
				*(dataIter+polIndex) = source[ch];
			}
			dataIter += polCount;
		}
	}
	else {
		switch(polSource) {
			case Polarization::StokesI: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsDest, polIndexA);
				bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsDest, polIndexB);
				if(!hasXX || !hasYY) {
					Polarization::TypeToIndex(Polarization::RR, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsDest, polIndexB);
				}
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					if(std::isfinite(source[ch].real()))
					{
						*(dataIter + polIndexA) = source[ch];
						*(dataIter + polIndexB) = source[ch];
					}
					dataIter += polCount;
				}
			}
			break;
			case Polarization::StokesQ: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXX = Polarization::TypeToIndex(Polarization::XX, polsDest, polIndexA);
				bool hasYY = Polarization::TypeToIndex(Polarization::YY, polsDest, polIndexB);
				if(hasXX && hasYY) {
					// StokesQ to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							casacore::Complex stokesI = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexB) + *(dataIter + polIndexA));
							*(dataIter + polIndexA) = stokesI - source[ch]; // XX = I - Q
							*(dataIter + polIndexB) = stokesI + source[ch]; // YY = I + Q
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesQ to circular
					Polarization::TypeToIndex(Polarization::RL, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							*(dataIter + polIndexA) = source[ch]; // rl = Q + iU (with U still zero)
							*(dataIter + polIndexB) = source[ch]; // lr = Q - iU (with U still zero)
						}
						dataIter += polCount;
					}
					throw std::runtime_error("Can't store polarization in set (not implemented or conversion not possible)");
				}
			}
			break;
			case Polarization::StokesU: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsDest, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsDest, polIndexB);
				if(hasXY && hasYX) {
					// StokesU to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							*(dataIter + polIndexA) = source[ch]; // XY = (U - iV), V still zero
							*(dataIter + polIndexB) = source[ch]; // YX = (U + iV), V still zero
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesU to circular
					Polarization::TypeToIndex(Polarization::RL, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// Q = (RL + LR) / 2
							casacore::Complex stokesQ = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexA) + *(dataIter + polIndexB));
							casacore::Complex iTimesStokesU = casacore::Complex(-source[ch].imag(), source[ch].real());
							*(dataIter + polIndexA) = stokesQ + iTimesStokesU; // rl = Q + iU
							*(dataIter + polIndexB) = stokesQ - iTimesStokesU; // lr = Q - iU
						}
						dataIter += polCount;
					}
					throw std::runtime_error("Can't store polarization in set (not implemented or conversion not possible)");
				}
			}
			break;
			case Polarization::StokesV: {
				size_t polIndexA=0, polIndexB=0;
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsDest, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsDest, polIndexB);
				if(hasXY && hasYX) {
					// StokesV to linear
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// U = (YX + XY)/2
							casacore::Complex stokesU = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexB) + *(dataIter + polIndexA));
							casacore::Complex iTimesStokesV = casacore::Complex(-source[ch].imag(), source[ch].real());
							*(dataIter + polIndexA) = stokesU - iTimesStokesV; // XY = (U - iV)
							*(dataIter + polIndexB) = stokesU + iTimesStokesV; // YX = (U + iV)
						}
						dataIter += polCount;
					}
				}
				else {
					// StokesV to circular
					Polarization::TypeToIndex(Polarization::RL, polsDest, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsDest, polIndexB);
					for(size_t ch=0; ch!=selectedChannelCount; ++ch)
					{
						if(std::isfinite(source[ch].real()))
						{
							// I = (RR + LL)/2
							casacore::Complex stokesI = casacore::Complex::value_type(0.5) * (*(dataIter + polIndexA) + *(dataIter + polIndexB));
							*(dataIter + polIndexA) = stokesI + source[ch]; // RR = I + V
							*(dataIter + polIndexB) = stokesI - source[ch]; // LL = I - V
						}
						dataIter += polCount;
					}
					throw std::runtime_error("Can't store polarization in set (not implemented or conversion not possible)");
				}
			}
			break;
			default:
				throw std::runtime_error("Can't store polarization in set (not implemented or conversion not possible)");
		}
	}
}

void MSProvider::getRowRange(casacore::MeasurementSet& ms, const MSSelection& selection, size_t& startRow, size_t& endRow)
{
	startRow = 0;
	endRow = ms.nrow();
	if(selection.HasInterval())
	{
		std::cout << "Determining first and last row index... " << std::flush;
		casacore::ROScalarColumn<double> timeColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::TIME));
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

void MSProvider::initializeModelColumn(casacore::MeasurementSet& ms)
{
	casacore::ROArrayColumn<casacore::Complex> dataColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::DATA));
	if(ms.isColumn(casacore::MSMainEnums::MODEL_DATA))
	{
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
		casacore::IPosition dataShape = dataColumn.shape(0);
		bool isDefined = modelColumn.isDefined(0);
		bool isSameShape = false;
		if(isDefined)
		{
			casacore::IPosition modelShape = modelColumn.shape(0);
			isSameShape = modelShape == dataShape;
		}
		if(!isDefined || !isSameShape)
		{
			std::cout << "WARNING: Your model column does not have the same shape as your data column: resetting MODEL column.\n";
			casacore::Array<casacore::Complex> zeroArray(dataShape);
			for(casacore::Array<casacore::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
				*i = std::complex<float>(0.0, 0.0);
			for(size_t row=0; row!=ms.nrow(); ++row)
				modelColumn.put(row, zeroArray);
		}
	}
	else { //if(!_ms.isColumn(casacore::MSMainEnums::MODEL_DATA))
		std::cout << "Adding model data column... " << std::flush;
		casacore::IPosition shape = dataColumn.shape(0);
		casacore::ArrayColumnDesc<casacore::Complex> modelColumnDesc(ms.columnName(casacore::MSMainEnums::MODEL_DATA), shape);
		try {
			ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
		} catch(std::exception& e)
		{
			ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
		}
		
		casacore::Array<casacore::Complex> zeroArray(shape);
		for(casacore::Array<casacore::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
			*i = std::complex<float>(0.0, 0.0);
		
		casacore::ArrayColumn<casacore::Complex> modelColumn(ms, casacore::MS::columnName(casacore::MSMainEnums::MODEL_DATA));
		
		for(size_t row=0; row!=ms.nrow(); ++row)
			modelColumn.put(row, zeroArray);
		
		std::cout << "DONE\n";
	}
}

vector<PolarizationEnum> MSProvider::GetMSPolarizations(casacore::MeasurementSet& ms)
{
	std::vector<PolarizationEnum> pols;
	casacore::MSPolarization polTable(ms.polarization());
	casacore::ROArrayColumn<int> corrTypeColumn(polTable, casacore::MSPolarization::columnName(casacore::MSPolarizationEnums::CORR_TYPE));
	casacore::Array<int> corrTypeVec(corrTypeColumn(0));
	for(casacore::Array<int>::const_contiter p=corrTypeVec.cbegin(); p!=corrTypeVec.cend(); ++p)
		pols.push_back(Polarization::AipsIndexToEnum(*p));
	return pols;
}
