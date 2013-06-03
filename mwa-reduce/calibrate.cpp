#include "calibrationmethod.h"
#include "model.h"
#include "banddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: calibrate <model> <measurementset.ms> <phases.txt> <gains.txt>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n"
			<< "Option -a will average over frequency before fitting, nr should specify the amount\n"
			<< "of desired channels.\n";
	} else {
		int argi = 1;
		if(argc <= argi + 2) throw std::runtime_error("Incorrect parameters");
		const char *modelName = argv[argi];
		const char *msName = argv[argi+1];
		const char *outNamePhases = argv[argi+2];
		const char *outNameGains = argv[argi+3];
		casa::MeasurementSet ms(msName);
		
		std::cout << "Reading model... " << std::flush;
		Model model(modelName);
		std::cout << "DONE\n";
		
		std::cout << "Reading meta data... " << std::flush;
		
		/**
		 * Read some meta data from the measurement set
		 */
		casa::MSAntenna aTable = ms.antenna();
		size_t antennaCount = aTable.nrow();
		
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		if(channelCount == 0) throw std::runtime_error("No channels in set");
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::MSField fieldTable = ms.field();
		casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		casa::Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
		casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
		casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
		
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		if(polarizationCount != 4)
			throw std::runtime_error("Pol count in MS != 4");
		
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
		
		std::cout << "DONE (" << timestepCount << ")\nReading data & predicting model... " << std::flush;
		
		CalibrationMethod calMethod(channelCount, antennaCount, timestepCount);
		Predicter predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount);
		predicter.Initialize(model);
		
		std::vector<std::complex<double> > modelValues(4 * channelCount);
		casa::Array<complex_t> data(dataShape);
		casa::Array<float> weights(dataShape);
		casa::Array<bool> flags(dataShape);
		size_t timeIndex = 0;
		time = timeColumn(0);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			if(timeColumn(rowIndex) != time)
			{
				++timeIndex;
				time = timeColumn(rowIndex);
				std::cout << '.' << std::flush;
			}
			// Cross correlation?
			size_t antenna1 = ant1Column.get(rowIndex), antenna2 = ant2Column.get(rowIndex);
			if(antenna1 != antenna2)
			{
				dataColumn.get(rowIndex, data);
				weightColumn.get(rowIndex, weights);
				flagColumn.get(rowIndex, flags);
				std::complex<float> *dataPtr = data.cbegin();
				float *weightsPtr = weights.cbegin();
				bool *flagPtr = flags.cbegin();
				
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					if(flagPtr[ch]) weightsPtr[ch] = 0.0;
					
					double lambda = bandData.ChannelWavelength(ch);
					std::complex<double> p = predicter.Predict(model, u/lambda, v/lambda, w/lambda, ch);
					modelValues[ch*4] = p;
					modelValues[ch*4 + 1] = 0;
					modelValues[ch*4 + 2] = 0;
					modelValues[ch*4 + 3] = p;
				}
					
				calMethod.AddData(dataPtr, weightsPtr, &modelValues[0], antenna1, antenna2, timeIndex);
			}
		}
		std::cout << "DONE\nCalibrating...\n";
		
		calMethod.Execute(0.0001);
	}
}
