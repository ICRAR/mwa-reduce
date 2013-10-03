#include "wsinversion.h"
#include "uvwdistribution.h"
#include "imageweights.h"
#include "buffered_lane.h"

#include <ms/MeasurementSets/MeasurementSet.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <tables/Tables/ArrColDesc.h>

#include <iostream>
#include <stdexcept>

#include <boost/thread/thread.hpp>

WSInversion::MSData::MSData() : matchingRows(0), totalRowsProcessed(0)
{ }

WSInversion::MSData::~MSData()
{ }

WSInversion::WSInversion() : InversionAlgorithm(), _hasFrequencies(false), _gridMode(LayeredImager::NearestNeighbour), _laneBufferSize(16)
{
	_cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
}
		
void WSInversion::initializeMeasurementSet(const string& measurementSet, WSInversion::MSData& msData)
{
	std::cout << "Opening " << measurementSet << "... " << std::flush;
	msData.ms.reset(new casa::MeasurementSet(measurementSet));
	casa::MeasurementSet& ms(*msData.ms);
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
	msData.bandData = BandData(ms.spectralWindow());
	if(HasChannelRange())
	{
		msData.startChannel = ChannelRangeStart();
		msData.endChannel = ChannelRangeEnd();
		if(msData.startChannel >= msData.bandData.ChannelCount() || msData.endChannel > msData.bandData.ChannelCount()
			|| msData.startChannel == msData.endChannel)
		{
			std::ostringstream str;
			str << "An invalid channel range was specified! Measurement set " << measurementSet << " only has " << msData.bandData.ChannelCount() << " channels, requested imaging range is " << msData.startChannel << " -- " << msData.endChannel << '.';
			throw std::runtime_error(str.str());
		}
	}
	else {
		msData.startChannel = 0;
		msData.endChannel = msData.bandData.ChannelCount();
	}
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	const BandData partBandData = BandData(msData.bandData, msData.startChannel, msData.endChannel);
	if(_hasFrequencies)
	{
		_freqLow = std::min(_freqLow, partBandData.LowestFrequency());
		_freqHigh = std::max(_freqHigh, partBandData.HighestFrequency());
		_bandStart = std::min(_bandStart, partBandData.BandStart());
		_bandEnd = std::max(_bandEnd, partBandData.BandEnd());
		_startTime = std::min(_startTime, timeColumn(0).getValue().get());
	} else {
		_freqLow = partBandData.LowestFrequency();
		_freqHigh = partBandData.HighestFrequency();
		_bandStart = partBandData.BandStart();
		_bandEnd = partBandData.BandEnd();
		_startTime = timeColumn(0).getValue().get();
		_hasFrequencies = true;
	}
	
	std::cout << 'C' << std::flush;
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	
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

	msData.rowStart = 0;
	msData.rowEnd = ms.nrow();
	if(HasInterval())
	{
		std::cout << "Determining first and last row index... " << std::flush;
		casa::MEpoch time = timeColumn(0);
		size_t timestepIndex = 0;
		for(size_t row = 0; row!=ms.nrow(); ++row)
		{
			if(time.getValue() != timeColumn(row).getValue())
			{
				++timestepIndex;
				if(timestepIndex == IntervalStart())
					msData.rowStart = row;
				if(timestepIndex == IntervalEnd())
				{
					msData.rowEnd = row;
					break;
				}
				time = timeColumn(row).getValue();
			}
		}
		std::cout << "DONE (" << msData.rowStart << '-' << msData.rowEnd << ")\n";
	}
	
	std::cout << "Determining min and max w... " << std::flush;
	msData.maxW= -1e100;
	msData.minW = 1e100;
	double maxBaseline = 0.0;
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			double uInM = uvwArray(0), vInM = uvwArray(1), wInM = uvwArray(2);
			double wHi = fabs(wInM / partBandData.SmallestWavelength());
			double wLo = fabs(wInM / partBandData.LongestWavelength());
			msData.maxW = std::max(msData.maxW, wHi);
			msData.minW = std::min(msData.minW, wLo);
			maxBaseline = std::max(maxBaseline, uInM*uInM + vInM*vInM + wInM*wInM);
		}
	}
	_beamSize = partBandData.SmallestWavelength() / sqrt(maxBaseline);
	std::cout << "DONE (w=[" << msData.minW << " -- " << msData.maxW << "] lambdas)\n";
}

void WSInversion::countSamplesPerLayer(MSData& msData)
{
	casa::MeasurementSet &ms(*msData.ms);
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	const BandData& bandData(msData.bandData);
	
	std::vector<size_t> sampleCount(WGridSize());
	msData.matchingRows = 0;
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			const double wInMeters = uvwArray(2);
			for(size_t ch=msData.startChannel; ch!=msData.endChannel; ++ch)
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
	std::unique_ptr<casa::ROArrayColumn<casa::Complex>> modelColumn;
	if(DoSubtractModel() && !DoImagePSF())
		modelColumn.reset(new casa::ROArrayColumn<std::complex<float>>(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
	
	casa::IPosition dataShape = dataColumn.shape(0);
	msData.polarizationCount = dataShape[0];
	
	const BandData selectedBand(msData.bandData, msData.startChannel, msData.endChannel);
	_imager->PrepareBand(selectedBand);
	
	casa::Array<std::complex<float>> msDataArr(dataShape), msModelDataArr(dataShape);
	casa::Array<float> msWeightsArr(dataShape);
	casa::Array<bool> msFlagsArr(dataShape);
	size_t rowsRead = 0;
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			const double
				wInMeters = uvwArray(2),
				w1 = wInMeters / selectedBand.LongestWavelength(),
				w2 = wInMeters / selectedBand.SmallestWavelength();
			if(_imager->IsInLayerRange(w1, w2))
			{
				InversionWorkItem newItem;
				newItem.u = uvwArray(0);
				newItem.v = uvwArray(1);
				newItem.w = wInMeters;
				newItem.data = new std::complex<float>[selectedBand.ChannelCount()];
				
				weightColumn.get(row, msWeightsArr);
				flagColumn.get(row, msFlagsArr);
				
				double rowWeight = 1.0;
				switch(Weighting())
				{
					case UniformWeighted:
					{
						float* msWeightIter = msWeightsArr.cbegin() + msData.startChannel * msData.polarizationCount;
						for(size_t ch=0; ch!=selectedBand.ChannelCount(); ++ch)
						{
							double
								u = newItem.u / selectedBand.ChannelWavelength(ch),
								v = newItem.v / selectedBand.ChannelWavelength(ch),
								weight = PrecalculatedWeightInfo()->GetUniformWeight(u, v);
							for(size_t p=0; p!=msData.polarizationCount; ++p)
							{
								*msWeightIter *= weight;
								++msWeightIter;
							}
						}
					} break;
					case NaturalWeighted:
						rowWeight = 1.0;
						break;
					case DistanceWeighted:
						rowWeight = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
						break;
				}
				
				dataColumn.get(row, msDataArr);
				if(DoImagePSF())
				{
					copyWeights(newItem.data, msData.startChannel, msData.endChannel, msData.polarizationCount, msDataArr, msWeightsArr, msFlagsArr, rowWeight);
				}
				else {
					if(DoSubtractModel())
					{
						modelColumn->get(row, msModelDataArr);
						casa::Array<casa::Complex>::contiter modelIter = msModelDataArr.cbegin();
						for(casa::Array<casa::Complex>::contiter iter = msDataArr.cbegin(); iter!=msDataArr.cend(); ++iter)
						{
							*iter -= *modelIter;
							modelIter++;
						}
					}
					copyWeightedData(newItem.data, msData.startChannel, msData.endChannel, msData.polarizationCount, msDataArr, msWeightsArr, msFlagsArr, rowWeight);
				}
				_inversionWorkLane->write(newItem);
				
				++rowsRead;
			}
		}
	}
	modelColumn.reset();
	std::cout << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsRead;
}

void WSInversion::workThreadParallel(const BandData* selectedBand)
{
	std::unique_ptr<lane<InversionWorkSample>[]> lanes(new lane<InversionWorkSample>[_cpuCount]);
	boost::thread_group group;
	// Samples of the same w-layer are collected in a buffer
	// before they are written into the lane. This is done because writing
	// to a lane is reasonably slow; it requires holding a mutex. Without
	// these buffers, writing the lane was a bottleneck and multithreading
	// did not help.
	std::unique_ptr<lane_write_buffer<InversionWorkSample>[]>
		bufferedLanes(new lane_write_buffer<InversionWorkSample>[_cpuCount]);
	for(size_t i=0; i!=_cpuCount; ++i)
	{
		lanes[i].resize(selectedBand->ChannelCount() * _laneBufferSize);
		bufferedLanes[i].reset(&lanes[i], std::max<size_t>(selectedBand->ChannelCount(), _laneBufferSize));
		
		group.add_thread(new boost::thread(&WSInversion::workThreadPerSample, this, &lanes[i]));
	}
	
	InversionWorkItem workItem;
	while(_inversionWorkLane->read(workItem))
	{
		InversionWorkSample sampleData;
		for(size_t ch=0; ch!=selectedBand->ChannelCount(); ++ch)
		{
			double wavelength = selectedBand->ChannelWavelength(ch);
			sampleData.sample = workItem.data[ch];
			sampleData.uInLambda = workItem.u / wavelength;
			sampleData.vInLambda = workItem.v / wavelength;
			sampleData.wInLambda = workItem.w / wavelength;
			size_t cpu = _imager->WToLayer(sampleData.wInLambda) % _cpuCount;
			bufferedLanes[cpu].write(sampleData);
		}
		delete[] workItem.data;
	}
	for(size_t i=0; i!=_cpuCount; ++i)
		bufferedLanes[i].write_end();
	group.join_all();
}

void WSInversion::workThreadPerSample(lane<InversionWorkSample>* workLane)
{
	lane_read_buffer<InversionWorkSample> buffer(workLane, std::min(_laneBufferSize*16, workLane->capacity()));
	InversionWorkSample sampleData;
	while(buffer.read(sampleData))
	{
		_imager->AddDataSample(sampleData.sample, sampleData.uInLambda, sampleData.vInLambda, sampleData.wInLambda);
	}
}


void WSInversion::sampleToMeasurementSet(MSData &msData)
{
	casa::MeasurementSet &ms(*msData.ms);
	ms.reopenRW();
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
	
	casa::ROArrayColumn<std::complex<float> > dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
	msData.polarizationCount = dataColumn.shape(0)[0];
	
	if(!ms.isColumn(casa::MSMainEnums::MODEL_DATA))
	{
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
		
		_modelColumn.reset(new casa::ArrayColumn<std::complex<float> >(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
		
		for(size_t row=0; row!=ms.nrow(); ++row)
			_modelColumn->put(row, zeroArray);
		
		std::cout << "DONE\n";
	}
	else {
		_modelColumn.reset(new casa::ArrayColumn<std::complex<float> >(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
	}
	
	const BandData selectedBandData(msData.bandData, msData.startChannel, msData.endChannel);
	_imager->PrepareBand(selectedBandData);
	
	size_t rowsProcessed = 0;
	
	lane<SamplingWorkItem> calcLane(_laneBufferSize+_cpuCount), writeLane(_laneBufferSize);
	boost::thread writeThread(&WSInversion::visSampleWriteThread, this, &writeLane, &msData);
	boost::thread_group calcThreads;
	for(size_t i=0; i!=_cpuCount; ++i)
		calcThreads.add_thread(new boost::thread(&WSInversion::visSampleCalcThread, this, &calcLane, &writeLane));

		
	/* Start by reading the u,v,ws in, so we don't need IO access
	 * from this thread during further processing */
	std::vector<double> us, vs, ws;
	std::vector<size_t> rowIndices;
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			casa::Vector<double> uvwArray = uvwColumn(row);
			const double
				wInMeters = uvwArray(2),
				w1 = wInMeters / selectedBandData.LongestWavelength(),
				w2 = wInMeters / selectedBandData.SmallestWavelength();
			if(_imager->IsInLayerRange(w1, w2))
			{
				us.push_back(uvwArray(0));
				vs.push_back(uvwArray(1));
				ws.push_back(uvwArray(2));
				rowIndices.push_back(row);
				++rowsProcessed;
			}
		}
	}
	
	for(size_t i=0; i!=us.size(); ++i)
	{
		SamplingWorkItem newItem;
		newItem.u = us[i];
		newItem.v = vs[i];
		newItem.w = ws[i];
		newItem.data = new std::complex<float>[selectedBandData.ChannelCount()];
		newItem.rowIndex = rowIndices[i];
				
		calcLane.write(newItem);
	}
	std::cout << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsProcessed;
	
	calcLane.write_end();
	calcThreads.join_all();
	writeLane.write_end();
	writeThread.join();
		
	_modelColumn.reset();
}

void WSInversion::visSampleCalcThread(lane<SamplingWorkItem>* inputLane, lane<SamplingWorkItem>* outputLane)
{
	SamplingWorkItem item;
	while(inputLane->read(item))
	{
		_imager->SampleData(item.data, item.u, item.v, item.w);
		
		outputLane->write(item);
	}
}

void WSInversion::visSampleWriteThread(lane<SamplingWorkItem>* samplingWorkLane, const MSData* msData)
{
	SamplingWorkItem workItem;
	
	// Read one, which makes sure other thread has initialized the ms
	if(!samplingWorkLane->read(workItem))
		return;
	
	casa::IPosition shape = _modelColumn->shape(0);
	casa::Array<std::complex<float>> data(shape);
	size_t polCount = shape[0];
	const BandData selectedBand = msData->SelectedBand();
	const size_t selectedChannelCount = selectedBand.ChannelCount();
	
	int polIndex = polarizationIndex();
	do
	{
		_modelColumn->get(workItem.rowIndex, data);
		casa::Array<std::complex<float>>::contiter dataIter = data.cbegin() + msData->startChannel * polCount;
		
		if(Polarization() == Polarization::StokesI)
		{
			for(size_t ch=0; ch!=selectedChannelCount; ++ch)
			{
				if(std::isfinite(workItem.data[ch].real()))
				{
					if(AddToModel())
					{
						*dataIter += workItem.data[ch];
						*(dataIter + (polCount-1)) += workItem.data[ch];
					} else {
						*dataIter = workItem.data[ch];
						*(dataIter + (polCount-1)) = workItem.data[ch];
					}
				}
				dataIter += polCount;
			}
		} else {
			for(size_t ch=0; ch!=selectedChannelCount; ++ch)
			{
				if(std::isfinite(workItem.data[ch].real()))
				{
					if(AddToModel())
						*(dataIter+polIndex) += workItem.data[ch];
					else
						*(dataIter+polIndex) = workItem.data[ch];
				}
				dataIter += polCount;
			}
		}
		_modelColumn->put(workItem.rowIndex, data);
	}
	while(samplingWorkLane->read(workItem));
}


void WSInversion::Invert()
{
	MSData msDataVector[MeasurementSetCount()];
	_hasFrequencies = false;
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

	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), _cpuCount));
	_imager->SetGridMode(_gridMode);
	_imager->PrepareWLayers(WGridSize(), memSize*2/4, minW, maxW);
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		countSamplesPerLayer(msDataVector[i]);
	
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		_inversionWorkLane.reset(new lane<InversionWorkItem>(16));
		
		_imager->StartInversionPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			MSData& msData = msDataVector[i];
			
			const BandData selectedBand(msData.SelectedBand());
			
			//boost::thread thread(&WSInversion::workThread, this, _inversionWorkLane.get());
			boost::thread thread(&WSInversion::workThreadParallel, this, &selectedBand);
		
			gridMeasurementSet(msData);
			
			_inversionWorkLane->write_end();
			thread.join();
		}
		
		std::cout << "Fourier transforming layers, w-term correction & summing in parallel...\n";
		_imager->FinishInversionPass();
	}
	
	size_t totalRowsRead = 0, totalMatchingRows = 0;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		totalRowsRead += msDataVector[i].totalRowsProcessed;
		totalMatchingRows += msDataVector[i].matchingRows;
	}
	
	std::cout << "Total rows read: " << totalRowsRead << " (overhead: " << round(totalRowsRead * 100.0 / totalMatchingRows - 100.0) << "%)\n";
	
	_imager->FinalizeImage(1.0/_totalWeight);
}

void WSInversion::InvertToVisibilities(const double *image)
{
	MSData msDataVector[MeasurementSetCount()];
	_hasFrequencies = false;
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

	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), _cpuCount));
	_imager->SetGridMode(_gridMode);
	_imager->PrepareWLayers(WGridSize(), memSize*2/4, minW, maxW);
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		countSamplesPerLayer(msDataVector[i]);
	
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "W-term correction & fourier transforming layers to uv space (in parallel)...\n";
		_imager->PrepareImageForVisibilitySampling(image, 1.0);
		
		std::cout << "Starting visibility sampling pass " << pass << ".\n";
		_imager->StartVisibilitySamplingPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			sampleToMeasurementSet(msDataVector[i]);
	}
	
	size_t totalRowsWritten = 0, totalMatchingRows = 0;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
	{
		totalRowsWritten += msDataVector[i].totalRowsProcessed;
		totalMatchingRows += msDataVector[i].matchingRows;
	}
	
	std::cout << "Total rows written: " << totalRowsWritten << " (overhead: " << round(totalRowsWritten * 100.0 / totalMatchingRows - 100.0) << "%)\n";
}

void WSInversion::copyWeightedData(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(Polarization() == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			weightPtr += polCount-1;
			inPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] += *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			++weightPtr;
			++inPtr;
			++flagPtr;
		}
	} else {
		int polIndex = polarizationIndex();
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			weightPtr += polCount;
			inPtr += polCount;
			flagPtr += polCount;
		}
	}
}

void WSInversion::copyWeights(std::complex<float>* dest, size_t startChannel, size_t endChannel, size_t polCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin() + startChannel * polCount;
	casa::Array<float>::const_contiter weightPtr = weights.cbegin() + startChannel * polCount;
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin() + startChannel * polCount;
	const size_t selectedChannelCount = endChannel - startChannel;
		
	if(Polarization() == Polarization::StokesI)
	{
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			inPtr += polCount-1;
			weightPtr += polCount-1;
			flagPtr += polCount-1;
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] += (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			++inPtr;
			++weightPtr;
			++flagPtr;
		}
	} else {
		int polIndex = polarizationIndex();
		
		inPtr += polIndex;
		weightPtr += polIndex;
		flagPtr += polIndex;
		for(size_t ch=0; ch!=selectedChannelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			inPtr += polCount;
			weightPtr += polCount;
			flagPtr += polCount;
		}
	}
}
