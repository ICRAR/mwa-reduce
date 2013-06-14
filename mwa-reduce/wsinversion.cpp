#include "wsinversion.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <iostream>
#include <stdexcept>

#include <boost/thread/thread.hpp>

void WSInversion::Execute()
{
	std::cout << "Opening " << MeasurementSetPath() << "... " << std::flush;
	casa::MeasurementSet ms(MeasurementSetPath());
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	/**
		* Read some meta data from the measurement set
		*/
	std::cout << 'A' << std::flush;
	casa::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	casa::MPosition ant1Pos = antPosColumn(0);
	
	std::cout << 'B' << std::flush;
	BandData bandData(ms.spectralWindow());
	size_t channelCount = bandData.ChannelCount();
	
	std::cout << 'C' << std::flush;
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<std::complex<float> > dataColumn(ms, DataColumnName());
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	
	std::cout << 'F' << std::flush;
	casa::MSField fTable(ms.field());
	if(fTable.nrow() != 1) throw std::runtime_error("Need exactly one field in set");
	casa::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	casa::MDirection refDir = refDirColumn(0);
	casa::MEpoch curtime = timeColumn(0);
	casa::MeasFrame frame(ant1Pos, curtime);
	casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MDirection j2000 = casa::MDirection::Convert(refDir, j2000Ref)();
	casa::Vector<casa::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];
	
	std::cout << 'D' << std::flush;
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	std::cout << " DONE (" << polarizationCount << ")\n";
	
	// Determine min and max w
	std::cout << "Determining min and max w... " << std::flush;
	double maxW= -1e100, minW = 1e100;
	for(size_t row=0;row!=ms.nrow();++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			double wInM = uvwArray(2);
			double wHi = fabs(wInM / bandData.SmallestWavelength());
			double wLo = fabs(wInM / bandData.LongestWavelength());
			maxW = std::max(maxW, wHi);
			minW = std::min(minW, wLo);
		}
	}
	std::cout << "DONE (min,max w=" << minW << ',' << maxW << " lambdas)\n";
	
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);
	std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";

	long cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), cpuCount));
	_imager->PrepareForObservation(WGridSize(), memSize*2/4, minW, maxW, bandData);
	
	std::vector<size_t> sampleCount(WGridSize());
	size_t matchingRows = 0;
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			const double wInMeters = uvwArray(2);
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				double w = wInMeters / bandData.ChannelWavelength(ch);
				++sampleCount[_imager->WToLayer(w)];
			}
			++matchingRows;
		}
	}
	std::cout << "Visibility count per layer: ";
	for(std::vector<size_t>::const_iterator i=sampleCount.begin(); i!=sampleCount.end(); ++i)
	{
		std::cout << *i << ' ';
	}
	std::cout << '\n';
	
	casa::Array<std::complex<float>> data(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	size_t totalRowsRead = 0;
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		_workLane.reset(new lane<WorkItem>(16));
		boost::thread thread(&WSInversion::workThread, this);
		
		_imager->StartPass(pass);
		
		size_t rowsRead = 0;
		for(size_t row=0; row!=ms.nrow(); ++row)
		{
			if(ant1Column(row) != ant2Column(row))
			{
				casa::Vector<double> uvwArray = uvwColumn(row);
				const double
					wInMeters = uvwArray(2),
					w1 = wInMeters / bandData.LongestWavelength(),
					w2 = wInMeters / bandData.SmallestWavelength();
				if(_imager->IsInLayerRange(w1, w2))
				{
					WorkItem newItem;
					newItem.u = uvwArray(0);
					newItem.v = uvwArray(1);
					newItem.w = wInMeters;
					newItem.data = new std::complex<float>[channelCount];
					
					weightColumn.get(row, weights);
					flagColumn.get(row, flags);
					
					if(DoImagePSF())
					{
						copyWeights(newItem.data, channelCount, weights, flags);
					}
					else {
						dataColumn.get(row, data);
						copyWeightedData(newItem.data, channelCount, data, weights, flags);
					}
					_workLane->write(newItem);
					
					++rowsRead;
				}
			}
		}
		std::cout << "Pass " << pass << ", rows that were required: " << rowsRead << '/' << matchingRows << '\n';
		totalRowsRead += rowsRead;
		
		_workLane->write_end();
		thread.join();
		
		std::cout << "Fourier transforming layers, w-term correction & summing in parallel...\n";
		_imager->FinishPass();
	}
	std::cout << "Total rows read: " << totalRowsRead << " (overhead: " << round(totalRowsRead * 100.0 / matchingRows - 100.0) << "%)\n";
	std::cout << "Total weight: " << _totalWeight << " Average per sample: " << _totalWeight / (totalRowsRead * channelCount) << '\n';
	_imager->FinalizeImage(1.0/_totalWeight);
}

void WSInversion::copyWeightedData(std::complex<float>* dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin();
	casa::Array<float>::const_contiter weightPtr = weights.cbegin();
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin();
		
	if(Polarization() == StokesI)
	{
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(*flagPtr)
				dest[ch] = 0;
			else {
				dest[ch] = *inPtr * *weightPtr;
				_totalWeight += *weightPtr;
			}
			weightPtr += 3;
			inPtr += 3;
			flagPtr += 3;
			if(!*flagPtr)
			{
				dest[ch] += *inPtr * *weightPtr;
				_totalWeight += *weightPtr;
			}
			++weightPtr;
			++inPtr;
			++flagPtr;
		}
	} else {
		int polIndex;
		switch(Polarization())
		{
			default: polIndex = 0; break;
			case XY: polIndex = 1; break;
			case YX: polIndex = 2; break;
			case YY: polIndex = 3; break;
		}
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(*flagPtr)
				dest[ch] = 0;
			else {
				dest[ch] = *inPtr * *weightPtr;
				_totalWeight += *weightPtr;
			}
			weightPtr += 4;
			inPtr += 4;
			flagPtr += 4;
		}
	}
}

void WSInversion::copyWeights(std::complex<float>* dest, size_t channelCount, const casa::Array<float>& weights, const casa::Array<bool>& flags)
{
	casa::Array<float>::const_contiter weightPtr = weights.cbegin();
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin();
		
	if(Polarization() == StokesI)
	{
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(*flagPtr)
				dest[ch] = 0;
			else {
				dest[ch] = *weightPtr;
				_totalWeight += *weightPtr;
			}
			weightPtr += 3;
			flagPtr += 3;
			if(!*flagPtr)
			{
				dest[ch] += *weightPtr;
				_totalWeight += *weightPtr;
			}
			++weightPtr;
			++flagPtr;
		}
	} else {
		int polIndex;
		switch(Polarization())
		{
			default: polIndex = 0; break;
			case XY: polIndex = 1; break;
			case YX: polIndex = 2; break;
			case YY: polIndex = 3; break;
		}
		
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(*flagPtr)
				dest[ch] = 0;
			else {
				dest[ch] = *weightPtr;
				_totalWeight += *weightPtr;
			}
			weightPtr += 4;
			flagPtr += 4;
		}
	}
}
