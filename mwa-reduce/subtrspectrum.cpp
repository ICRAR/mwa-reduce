#include "model.h"
#include "banddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <iostream>

int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		std::cout << "Syntax: subtrspectrum <model> <ms>\nDetermines spectrum for each source in the model at each time step,\nand subtracts it from the measurement set.";
		return -1;
	}
	
	Model model(argv[1]);
	casa::MeasurementSet ms(argv[2]);
	BandData bandData(ms.spectralWindow());
	
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Need 4 polarizations");
	
	casa::Array<casa::Complex> dataArray(dataShape);
	casa::Array<bool> flagArray(dataShape);
	casa::Array<float> weightArray(dataShape);
	
	double time = timeColumn(0);
	size_t timeStepStartRow = 0;
	
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		double thisTime = timeColumn(row);
		if(thisTime != time)
		{
			// Subtract data
			
			time = thisTime;
			timeStepStartRow = row+1;
		}
		
		
		dataColumn.get(row, dataArray);
		flagColumn.get(row, flagArray);
		weightColumn.get(row, weightArray);
		
		// Push data per source into thread
	}
}
