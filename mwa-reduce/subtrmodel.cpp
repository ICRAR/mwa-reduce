#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "banddata.h"
#include "sourcesdf.h"
#include "model.h"
#include "predicter.h"

using namespace casa;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: subtrmodel [-r] <model> <ms>\n"
			"Subtracts the model from the visibilities. This 'peels' the\n"
			"sources out. Only affects cross-correlations. -r to revert or -s to set.\n";
	} else {
		bool revert = false , setvis = false;
		size_t argi = 1;
		if(strcmp(argv[1], "-r") == 0) { revert=true; ++argi; }
		else if(strcmp(argv[1], "-s") == 0) { setvis=true; ++argi; }
		std::cout << "Reading model... " << std::flush;
		Model model(argv[argi]);
		std::cout << "DONE\n";
		
		std::cout << "Opening measurement set... " << std::flush;
		MeasurementSet ms(argv[argi+1], Table::Update);
		
		/**
		 * Read some meta data from the measurement set
		 */
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		
		MSField fieldTable = ms.field();
		ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
		ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
		ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		ROArrayColumn<double> uvwColumn(ms, ms.columnName(MSMainEnums::UVW));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		std::cout << "DONE\n";
		std::cout << "RA=" << phaseCentreRA << " Dec=" << phaseCentreDec << '\n';
		
		Predicter predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount);
		predicter.Initialize(model);
		
		/**
		 * Subtract
		 */
		if(revert)
			std::cout << "Adding back ";
		else
			std::cout << "Subtracting ";
		std::cout << model.SourceCount() << " sources... " << std::flush;
		Array<complex_t> data(dataShape);
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				dataColumn.get(rowIndex, data);
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				
				Array<complex_t>::iterator dataPtr = data.begin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					double lambda = bandData.ChannelWavelength(ch);
					if(setvis)
					{
						for(size_t p=0;p!=polarizationCount;++p) {
							*dataPtr = predicter.Predict(model, u/lambda, v/lambda, w/lambda, ch, p);
							++dataPtr;
						}
					} else {
						for(size_t p=0;p!=polarizationCount;++p)
						{
							std::complex<float> predicted;
							if(revert)
								predicted = predicter.Predict(model, u/lambda, v/lambda, w/lambda, ch, p);
							else
								predicted = -predicter.Predict(model, u/lambda, v/lambda, w/lambda, ch, p);
							*dataPtr += predicted;
							++dataPtr;
						}
					}
				}
				dataColumn.put(rowIndex, data);
				//std::cout << '.' << std::flush;
			}
		}
		
		std::cout << "DONE\n";		
	}
}
