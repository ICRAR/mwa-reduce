#include "msprovider.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrColDesc.h>

#include "../msselection.h"

void MSProvider::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
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
					throw std::runtime_error("can not form requested polarization (Stokes I) from available polarizations");
			}
			
			for(size_t ch=0; ch!=selectedChannelCount; ++ch)
			{
				weightPtr += polIndexA;
				inPtr += polIndexA;
				flagPtr += polIndexA;
				
				bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
				casa::Complex valA = *inPtr * (*weightPtr);
				
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
					casa::Complex valA = *inPtr * (*weightPtr);
					
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
					throw std::runtime_error("can not form requested polarization (Stokes Q) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casa::Complex valA = *inPtr * (*weightPtr);
					
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
				// Convert to StokesU from XX and YY
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casa::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else
						dest[ch] = valA + (*inPtr * (*weightPtr));
					
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
					throw std::runtime_error("can not form requested polarization (Stokes U) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casa::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casa::Complex diff = (valA - *inPtr * (*weightPtr));
						// U = -i (RL - LR)
						dest[ch] = casa::Complex(diff.imag(), -diff.real());
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
					casa::Complex valA = *inPtr * (*weightPtr);
					
					weightPtr += polIndexB - polIndexA;
					inPtr += polIndexB - polIndexA;
					flagPtr += polIndexB - polIndexA;
					
					bool flagB = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					if(flagA || flagB)
						dest[ch] = 0.0;
					else {
						casa::Complex diff = (valA - *inPtr * (*weightPtr));
						// V = -i (YX - XY)
						dest[ch] = casa::Complex(diff.imag(), -diff.real());
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
					throw std::runtime_error("can not form requested polarization (Stokes U) from available polarizations");
				for(size_t ch=0; ch!=selectedChannelCount; ++ch)
				{
					weightPtr += polIndexA;
					inPtr += polIndexA;
					flagPtr += polIndexA;
					
					bool flagA = *flagPtr || !std::isfinite(inPtr->real())|| !std::isfinite(inPtr->imag());
					casa::Complex valA = *inPtr * (*weightPtr);
					
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
void MSProvider::copyWeights(NumType* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut)
{
	const size_t polCount = polsIn.size();
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
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
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexA);
				}
			}
			break;
			case Polarization::StokesQ: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XX, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YY, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexA);
				}
			}
			break;
			case Polarization::StokesU: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RL, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LR, polsIn, polIndexA);
				}
			}
			break;
			case Polarization::StokesV: {
				bool hasXY = Polarization::TypeToIndex(Polarization::XY, polsIn, polIndexA);
				bool hasYX = Polarization::TypeToIndex(Polarization::YX, polsIn, polIndexB);
				if(!hasXY || !hasYX) {
					Polarization::TypeToIndex(Polarization::RR, polsIn, polIndexA);
					Polarization::TypeToIndex(Polarization::LL, polsIn, polIndexA);
				}
			}
			break;
			default:
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
void MSProvider::copyWeights<float>(float* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

template
void MSProvider::copyWeights<std::complex<float>>(std::complex<float>* dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum>& polsIn, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, PolarizationEnum polOut);

void MSProvider::reverseCopyData(casa::Array<std::complex<float>>& dest, size_t startChannel, size_t endChannel, const std::vector<PolarizationEnum> &polsDest, const std::complex<float>* source, PolarizationEnum polSource)
{
	size_t polCount = polsDest.size();
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

void MSProvider::initializeModelColumn(casa::MeasurementSet& ms)
{
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, casa::MS::columnName(casa::MSMainEnums::DATA));
	if(ms.isColumn(casa::MSMainEnums::MODEL_DATA))
	{
		casa::ArrayColumn<casa::Complex> modelColumn(ms, casa::MS::columnName(casa::MSMainEnums::MODEL_DATA));
		casa::IPosition dataShape = dataColumn.shape(0);
		bool isDefined = modelColumn.isDefined(0);
		bool isSameShape = false;
		if(isDefined)
		{
			casa::IPosition modelShape = modelColumn.shape(0);
			isSameShape = modelShape == dataShape;
		}
		if(!isDefined || !isSameShape)
		{
			std::cout << "WARNING: Your model column does not have the same shape as your data column: resetting MODEL column.\n";
			casa::Array<casa::Complex> zeroArray(dataShape);
			for(casa::Array<casa::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
				*i = std::complex<float>(0.0, 0.0);
			for(size_t row=0; row!=ms.nrow(); ++row)
				modelColumn.put(row, zeroArray);
		}
	}
	else { //if(!_ms.isColumn(casa::MSMainEnums::MODEL_DATA))
		std::cout << "Adding model data column... " << std::flush;
		casa::IPosition shape = dataColumn.shape(0);
		casa::ArrayColumnDesc<casa::Complex> modelColumnDesc(ms.columnName(casa::MSMainEnums::MODEL_DATA), shape);
		try {
			ms.addColumn(modelColumnDesc, "StandardStMan", true, true);
		} catch(std::exception& e)
		{
			ms.addColumn(modelColumnDesc, "StandardStMan", false, true);
		}
		
		casa::Array<casa::Complex> zeroArray(shape);
		for(casa::Array<casa::Complex>::contiter i=zeroArray.cbegin(); i!=zeroArray.cend(); ++i)
			*i = std::complex<float>(0.0, 0.0);
		
		casa::ArrayColumn<casa::Complex> modelColumn(ms, casa::MS::columnName(casa::MSMainEnums::MODEL_DATA));
		
		for(size_t row=0; row!=ms.nrow(); ++row)
			modelColumn.put(row, zeroArray);
		
		std::cout << "DONE\n";
	}
}

vector<PolarizationEnum> MSProvider::GetMSPolarizations(casa::MeasurementSet& ms)
{
	std::vector<PolarizationEnum> pols;
	casa::MSPolarization polTable(ms.polarization());
	casa::ROArrayColumn<int> corrTypeColumn(polTable, casa::MSPolarization::columnName(casa::MSPolarizationEnums::CORR_TYPE));
	casa::Array<int> corrTypeVec(corrTypeColumn(0));
	for(casa::Array<int>::const_contiter p=corrTypeVec.cbegin(); p!=corrTypeVec.cend(); ++p)
		pols.push_back(Polarization::AipsIndexToEnum(*p));
	return pols;
}
