#include "calibrationmethod.h"
#include "model.h"
#include "banddata.h"
#include "solutionfile.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <memory>

int main(int argc, char *argv[])
{
	if(argc < 4)
	{
		std::cout
			<< "Usage: calibrate [-p <phases.txt> <gains.txt>] [-l <limit>] [-i <niter>] [-m <model>] <measurementset.ms> <solutions.bin>\n\n"
			<< "This will calculate \"static\" phase offsets for all stations. It produces approximate least-squares solutions.\n"
			<< "Option -a will average over frequency before fitting, nr should specify the amount\n"
			<< "of desired channels.\n";
	} else {
		int argi = 1;
		bool savePlotFiles = false;
		std::string plotPhaseFile, plotGainFile, modelFile;
		size_t niter = 25;
		double limit = 0.0001;
		
		while(argv[argi][0] == '-')
		{
			if(strcmp(argv[argi], "-p") == 0)
			{
				savePlotFiles = true;
				plotPhaseFile = argv[argi+1];
				plotGainFile = argv[argi+2];
				argi += 3;
			}
			else if(strcmp(argv[argi], "-i") == 0)
			{
				niter = atoi(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-l") == 0)
			{
				limit = atof(argv[argi+1]);
				argi += 2;
			}
			else if(strcmp(argv[argi], "-m") == 0)
			{
				modelFile = argv[argi+1];
				argi += 2;
			}
			else throw std::runtime_error("Invalid parameter");
		}
		
		if(argc <= argi + 1) throw std::runtime_error("Incorrect parameters");
		
		const char *msName = argv[argi];
		const char *outName = argv[argi+1];
		casa::MeasurementSet ms(msName);
		
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
		std::unique_ptr<casa::ROArrayColumn<complex_t>> modelColumn;
		
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
		
		CalibrationMethod calMethod(channelCount, antennaCount, timestepCount);
		std::cout << "DONE (" << timestepCount << ")\n";
		
		std::unique_ptr<Model> model;
		std::unique_ptr<Predicter> predicter;
		if(modelFile.empty())
		{
			std::cout << "Reading data and model column... " << std::flush;
			modelColumn.reset(new casa::ROArrayColumn<complex_t>(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
		}
		else {
			std::cout << "Reading model... " << std::flush;
			model.reset(new Model(modelFile.c_str()));
			std::cout << "DONE\n";
			predicter.reset(new Predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount));
			predicter->Initialize(*model);
			predicter->ReportSources(*model);
			std::cout << "Reading data & predicting model... " << std::flush;
		}
		
		std::vector<std::complex<double> > modelValues(4 * channelCount);
		casa::Array<complex_t> data(dataShape), modelData(dataShape);
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
				
				if(modelFile.empty())
				{
					modelColumn->get(rowIndex, modelData);
					std::complex<float> *modelDataPtr = modelData.cbegin();
					for(size_t ch = 0; ch!=channelCount; ++ch)
					{
						for(size_t p=0; p!=4; ++p)
						{
							modelValues[ch*4+p] = *modelDataPtr;
							++modelDataPtr;
						}
					}
				}
				else {
					for(size_t ch = 0; ch!=channelCount; ++ch)
					{
						double lambda = bandData.ChannelWavelength(ch);
						for(size_t p=0; p!=4; ++p)
						{
							if(flagPtr[ch*4+p]) weightsPtr[ch*4+p] = 0.0;
							std::complex<double> pVal = predicter->Predict(*model, u/lambda, v/lambda, w/lambda, ch, p);
							modelValues[ch*4+p] = pVal;
						}
					}
				}
					
				calMethod.AddData(dataPtr, weightsPtr, &modelValues[0], antenna1, antenna2, timeIndex);
			}
		}
		std::cout << "DONE\nCalibrating...\n";
		
		calMethod.Execute(limit, niter);

		SolutionFile solutionFile;
		solutionFile.SetAntennaCount(antennaCount);
		solutionFile.SetChannelCount(channelCount);
		solutionFile.SetPolarizationCount(4);
		solutionFile.OpenForWriting(outName);

		for(size_t ant=0; ant!=antennaCount; ++ant)
		{
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				for(size_t p=0; p!=4; ++p)
				{
					const std::complex<double> val = calMethod.JonesSolution(ant, ch, p);
					solutionFile.WriteNextSolution(val);
				}
			}
		}
		
		
		if(savePlotFiles)
		{
			std::ofstream phasePlotStream(plotPhaseFile.c_str()), gainPlotStream(plotGainFile.c_str());
			phasePlotStream << antennaCount << ' ' << channelCount << " 4\n";
			gainPlotStream << antennaCount << ' ' << channelCount << " 4\n";
			
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				phasePlotStream << ch << '\t';
				gainPlotStream << ch << '\t';
				for(size_t p=0; p!=4; ++p)
				{
					for(size_t ant=0; ant!=antennaCount; ++ant)
					{
						const std::complex<double> val = calMethod.JonesSolution(ant, ch, p);
						double s1, s2;
						calMethod.SolutionSingularValue(ant, ch, s1, s2);
						switch(p)
						{
							case 0: gainPlotStream << '\t' << s1; break;
							case 1: case 2: gainPlotStream << '\t' << 0.0; break;
							case 3: gainPlotStream << '\t' << s2;
						}
						phasePlotStream << '\t' << std::arg(val);
					}
				}
				phasePlotStream << '\n';
				gainPlotStream << '\n';
			}
		}
	}
}
