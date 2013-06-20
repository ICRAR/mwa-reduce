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

void WSInversion::initializeMeasurementSet(const string& measurementSet, WSInversion::MSData& msData)
{
	std::cout << "Opening " << measurementSet << "... " << std::flush;
	msData.ms = new casa::MeasurementSet(measurementSet);
	casa::MeasurementSet &ms(*msData.ms);
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
	msData.channelCount = bandData.ChannelCount();
	_freqHigh = bandData.HighestFrequency();
	_freqLow = bandData.LowestFrequency();
	
	std::cout << 'C' << std::flush;
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
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
	std::cout << " DONE\n";
	
	// Determine min and max w
	std::cout << "Determining min and max w... " << std::flush;
	msData.maxW= -1e100;
	msData.minW = 1e100;
	double maxBaseline = 0.0;
	for(size_t row=0;row!=ms.nrow();++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			double uInM = uvwArray(0), vInM = uvwArray(1), wInM = uvwArray(2);
			double wHi = fabs(wInM / bandData.SmallestWavelength());
			double wLo = fabs(wInM / bandData.LongestWavelength());
			msData.maxW = std::max(msData.maxW, wHi);
			msData.minW = std::min(msData.minW, wLo);
			maxBaseline = std::max(maxBaseline, uInM*uInM + vInM*vInM + wInM*wInM);
		}
	}
	_beamSize = bandData.SmallestWavelength() / sqrt(maxBaseline);
	std::cout << "DONE (min,max w=" << msData.minW << ',' << msData.maxW << " lambdas)\n";
}

void WSInversion::countSamplesPerLayer(MSData& msData)
{
	casa::MeasurementSet &ms(*msData.ms);
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	BandData bandData(ms.spectralWindow());
	
	std::vector<size_t> sampleCount(WGridSize());
	msData.matchingRows = 0;
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			const double wInMeters = uvwArray(2);
			for(size_t ch=0; ch!=msData.channelCount; ++ch)
			{
				double w = wInMeters / bandData.ChannelWavelength(ch);
				++sampleCount[_imager->WToLayer(w)];
			}
			++msData.matchingRows;
		}
	}
	std::cout << "Visibility count per layer: ";
	for(std::vector<size_t>::const_iterator i=sampleCount.begin(); i!=sampleCount.end(); ++i)
	{
		std::cout << *i << ' ';
	}
	std::cout << '\n';
}

void WSInversion::gridMeasurementSet(MSData &msData)
{
	casa::MeasurementSet &ms(*msData.ms);
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<std::complex<float> > dataColumn(ms, DataColumnName());
	
	casa::IPosition dataShape = dataColumn.shape(0);
	msData.polarizationCount = dataShape[0];
	
	BandData bandData(ms.spectralWindow());
	_imager->PrepareBand(bandData);
	
	casa::Array<std::complex<float>> data(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
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
				newItem.data = new std::complex<float>[msData.channelCount];
				
				weightColumn.get(row, weights);
				flagColumn.get(row, flags);
				
				double rowWeight;
				switch(Weighting())
				{
					default:
					case NaturalWeighted:
						rowWeight = 1.0;
						break;
					case DistanceWeighted:
						rowWeight = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
						break;
				}
				
				if(DoImagePSF())
				{
					copyWeights(newItem.data, msData.channelCount, weights, flags, rowWeight);
				}
				else {
					dataColumn.get(row, data);
					copyWeightedData(newItem.data, msData.channelCount, data, weights, flags, rowWeight);
				}
				_workLane->write(newItem);
				
				++rowsRead;
			}
		}
	}
	std::cout << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsRead += rowsRead;
}

void WSInversion::Execute()
{
	std::vector<MSData> msDataVector(MeasurementSetCount());
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		initializeMeasurementSet(MeasurementSetPath(i), msDataVector[i]);
	
	double minW = msDataVector[0].minW;
	double maxW = msDataVector[0].maxW;
	for(size_t i=1; i!=MeasurementSetCount(); ++i)
	{
		if(msDataVector[i].minW < minW) minW = msDataVector[i].minW;
		if(msDataVector[i].maxW > maxW) maxW = msDataVector[i].maxW;
	}
	
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	int64_t memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) memSize / (1024.0*1024.0*1024.0);
	std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";

	long cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), cpuCount));
	_imager->PrepareWLayers(WGridSize(), memSize*2/4, minW, maxW);
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		countSamplesPerLayer(msDataVector[i]);
	
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		_workLane.reset(new lane<WorkItem>(16));
		boost::thread thread(&WSInversion::workThread, this);
		
		_imager->StartPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			gridMeasurementSet(msDataVector[i]);
		
		_workLane->write_end();
		thread.join();
		
		std::cout << "Fourier transforming layers, w-term correction & summing in parallel...\n";
		_imager->FinishPass();
	}
	
	size_t totalRowsRead = 0, totalMatchingRows = 0;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		totalRowsRead += msDataVector[i].totalRowsRead;
		totalMatchingRows += msDataVector[i].matchingRows;
	}
	
	std::cout << "Total rows read: " << totalRowsRead << " (overhead: " << round(totalRowsRead * 100.0 / totalMatchingRows - 100.0) << "%)\n";
	
	_imager->FinalizeImage(1.0/_totalWeight);
}

void WSInversion::copyWeightedData(std::complex<float>* dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
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
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 3;
			inPtr += 3;
			flagPtr += 3;
			if(!*flagPtr)
			{
				dest[ch] += *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
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
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 4;
			inPtr += 4;
			flagPtr += 4;
		}
	}
}

void WSInversion::copyWeights(std::complex<float>* dest, size_t channelCount, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
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
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 3;
			flagPtr += 3;
			if(!*flagPtr)
			{
				dest[ch] += (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
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
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			weightPtr += 4;
			flagPtr += 4;
		}
	}
}
