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
		std::cout << "Usage: spectrum <model> <ms>\n"
			"Calculates the spectrum directly from the ms, for each source in the model.\n";
	} else {
		size_t argi = 1;
		Model model(argv[argi]);
		Model measuredModel(model);
		// Set all sources to flux 1 Jy
		for(Model::iterator source=model.begin();source!=model.end();++source)
		{
			ModelSource &s = *source;
			s.SetBrightness(SourceSDFWithSI<long double>(1.0, 0.0, 1.0));
		}
		
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
		ArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
		ROArrayColumn<double> uvwColumn(ms, ms.columnName(MSMainEnums::UVW));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		Predicter predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount);
		predicter.Initialize(model);
		
		/**
		 * Calculate spectra
		 */
		Array<complex_t> data(dataShape);
		Array<bool> flags(dataShape);
		long double sourceFlux[model.SourceCount()*channelCount];
		long unsigned sourceMeasCount[model.SourceCount()*channelCount];
		for(size_t i=0;i!=model.SourceCount()*channelCount;++i)
		{
			sourceFlux[i] = 0.0;
			sourceMeasCount[i] = 0;
		}
		
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				
				Array<complex_t>::iterator dataPtr = data.begin();
				Array<bool>::iterator flagPtr = flags.begin();
				for(size_t ch=0; ch!=channelCount; ++ch)
				{
					double lambda = bandData.ChannelWavelength(ch);
					for(size_t p=0;p!=polarizationCount;++p)
					{
						float real = dataPtr->real(), imag = dataPtr->imag();
						if(!(*flagPtr) && std::isfinite(real) && std::isfinite(imag) && (polarizationCount!=4 || p==0 || p==3))
						{
							size_t sourceIndex = ch * model.SourceCount();
							for(Model::const_iterator source=model.begin();source!=model.end();++source)
							{
								Predicter::CNumType predicted = predicter.Predict(*source, u/lambda, v/lambda, w/lambda, ch);
								// add real(data * conj(predicted))
								sourceFlux[sourceIndex] += real * predicted.real() + imag * predicted.imag();
								sourceMeasCount[sourceIndex]++;
								++sourceIndex;
							}
						}
						++dataPtr;
						++flagPtr;
					}
				}
			}
		}
		
		bool outputModel = true;
		if(outputModel)
		{
			size_t sourceIndex = 0;
			for(Model::iterator source=model.begin();source!=model.end();++source)
			{
				SourceSDFWithSamples<long double> sdf;
				size_t itemIndex = sourceIndex;
				for(size_t ch=0; ch!=channelCount;++ch)
				{
					sdf.AddSample((sourceFlux[itemIndex] / sourceMeasCount[itemIndex]), bandData.ChannelFrequency(ch));
					itemIndex += model.SourceCount();
				}
				source->SetBrightness(sdf);
				std::cout << source->ToStringLine() << '\n';
				++sourceIndex;
			}
		} else {
			float sums[model.SourceCount()];
			size_t counts[model.SourceCount()];
			for(size_t i=0;i!=model.SourceCount();++i)
			{
				sums[i] = 0;
				counts[i] = 0;
			}
			
			for(size_t ch=0; ch!=channelCount;++ch)
			{
				std::cout << ch << '\t' << bandData.ChannelFrequency(ch);
				size_t sourceIndex = ch * model.SourceCount();
				for(size_t s=0; s!=model.SourceCount();++s)
				{
					std::cout << '\t' << (sourceFlux[sourceIndex] / sourceMeasCount[sourceIndex]);
					sums[s] += sourceFlux[sourceIndex];
					counts[s] += sourceMeasCount[sourceIndex];
					++sourceIndex;
				}
				std::cout << '\n';
			}
			std::cout << "avg\tavg\t";
			for(size_t i=0;i!=model.SourceCount();++i)
			{
				std::cout << '\t' << (sums[i] / counts[i]);
			}
			std::cout << '\n';
		}
	}
}
