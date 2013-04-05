#include "gausencoder.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <iostream>
#include <stdexcept>

typedef long double num_t;
typedef std::complex<float> complex_t;

using namespace casa;

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cout << "Usage: simcompress <quantcount> <stddev> <ms>\n";
		return -1;
	}
	
	size_t quantCount = atoi(argv[1]);
	num_t stddev = atof(argv[2]);
	const char *filename = argv[3];
	
	MeasurementSet ms(filename, Table::Update);
	
	/**
		* Read some meta data from the measurement set
		*/
	std::cout << "Opening measurement set... " << std::flush;
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
	ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
	ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
	ArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
	
	IPosition dataShape = dataColumn.shape(0);
	std::cout << "DONE\n";
	
	std::cout << "Initializing encoder... " << std::flush;
	GausEncoder<float> encoder(quantCount, stddev);
	std::cout << "DONE\n" << std::flush;
	
	/**
		* Apply compression effects
		*/
	std::cout << "Applying compression... " << std::flush;
	Array<complex_t> data(dataShape);
	Array<bool> flags(dataShape);
	num_t realError = 0.0, imagError = 0.0;
	num_t flaggedRealError = 0.0, flaggedImagError = 0.0;
	size_t count = 0, flaggedCount = 0;
	for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
	{
		//size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
		dataColumn.get(rowIndex, data);
		flagColumn.get(rowIndex, flags);
		Array<complex_t>::iterator dataPtr = data.begin();
		Array<bool>::const_iterator flagPtr = flags.begin();
		
		while(dataPtr != data.end())
		{
			if(std::isfinite(dataPtr->real()) && std::isfinite(dataPtr->imag()))
			{
				float newReal = encoder.Decode(encoder.Encode(dataPtr->real()));
				float newImag = encoder.Decode(encoder.Encode(dataPtr->imag()));
				
				num_t errReal = (newReal - dataPtr->real()) * (newReal - dataPtr->real());
				num_t errImag = (newImag - dataPtr->imag()) * (newImag - dataPtr->imag());
				
				if(*flagPtr)
				{
					flaggedRealError += errReal;
					flaggedImagError += errImag;
					++flaggedCount;
				} else {
					realError += errReal;
					imagError += errImag;
					++count;
				}
				
				*dataPtr = complex_t(newReal, newImag);
			}
			
			++dataPtr;
			++flagPtr;
		}
		
		dataColumn.put(rowIndex, data);
	}
	
	std::cout << "DONE\n";
	
	std::cout << "Statics on non-flagged data:\n";
	std::cout << "RMS on reals: " << sqrtl(realError / count) << "  RMS on imags: " << sqrtl(imagError / count) << '\n';
	
	std::cout << "Statics on flagged data:\n";
	std::cout << "RMS on reals: " << sqrtl(flaggedRealError / flaggedCount) << "  RMS on imags: " << sqrtl(flaggedImagError / flaggedCount) << '\n';
	
	return 0;
}
