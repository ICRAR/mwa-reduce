#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "beamevaluator.h"
#include "banddata.h"
#include "sourcesdf.h"
#include "model.h"
#include "predicter.h"
#include "mspredicter.h"
#include "progressbar.h"

using namespace casa;

template<typename T>
void addGausNoise(std::complex<T> &value, double sigma)
{
	long double x1, x2, w;

	do {
		long double r1 = (long double) rand() / (long double) RAND_MAX; 
		long double r2 = (long double) rand() / (long double) RAND_MAX; 
		x1 = 2.0 * r1 - 1.0;
		x2 = 2.0 * r2 - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = std::sqrt( (-2.0 * std::log( w ) ) / w ) * sigma;
	value += std::complex<T>(x1 * w, x2 * w);
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: subtrmodel [-r] <model> <ms>\n"
			"Subtracts the model from the visibilities. This 'peels' the\n"
			"sources out. Only affects cross-correlations. -r to revert or -s to set.\n";
	} else {
		bool revert = false , setvis = false, addNoise = false;
		double noiseSigma = 1.0;
		size_t argi = 1;
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-r") == 0) { revert=true; }
			else if(strcmp(argv[argi], "-s") == 0) { setvis=true; }
			else if(strcmp(argv[argi], "-n") == 0) { addNoise=true; ++argi; noiseSigma = atof(argv[argi]); }
			else throw std::runtime_error("Invalid param");
			++argi;
		}
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
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		ArrayColumn<complex_t> dataColumn(ms, ms.columnName(MSMainEnums::DATA));
		
		IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		std::cout << "DONE\n";
		
		BeamEvaluator beamEvaluator(ms);
		
		MSPredicter predicter(ms, model);
		predicter.Start(true);
		
		/**
		 * Subtract
		 */
		std::ostringstream taskDesc;
		if(revert)
			taskDesc << "Adding back ";
		else if(setvis)
			taskDesc << "Setting visibilities from ";
		else
			taskDesc << "Subtracting ";
		taskDesc << model.SourceCount() << " sources... " << std::flush;
		ProgressBar progress(taskDesc.str());
		
		Array<complex_t> data(dataShape);
		MSPredicter::RowData rowData;
		while(predicter.GetNextRow(rowData))
		{
			size_t rowIndex = rowData.rowIndex;
			
			boost::mutex::scoped_lock lock(predicter.IOMutex());
			progress.SetProgress(rowIndex, ms.nrow());
			dataColumn.get(rowIndex, data);
			lock.unlock();
			
			Array<complex_t>::iterator dataPtr = data.begin();
			std::complex<double> *modelDataPtr = rowData.modelData;
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					std::complex<double> predicted;
					if(revert || setvis)
						predicted = *modelDataPtr;
					else
						predicted = -*modelDataPtr;
					if(addNoise)
						addGausNoise(predicted, noiseSigma);
					if(setvis)
						*dataPtr = predicted;
					else
						*dataPtr += predicted;
					++dataPtr;
					++modelDataPtr;
				}
			}
			
			lock.lock();
			dataColumn.put(rowIndex, data);
			lock.unlock();
			
			predicter.FinishRow(rowData);
		}
		
		std::cout << "DONE\n";
	}
}
