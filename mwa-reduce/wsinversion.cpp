#include "wsinversion.h"
#include "uvwdistribution.h"
#include "imageweights.h"
#include "buffered_lane.h"
#include "partitionedms.h"

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

WSInversion::WSInversion(ImageBufferAllocator<double>* imageAllocator) : InversionAlgorithm(), _phaseCentreDL(0.0), _phaseCentreDM(0.0), _denormalPhaseCentre(false), _hasFrequencies(false), _gridMode(LayeredImager::NearestNeighbour), _laneBufferSize(16), _imageBufferAllocator(imageAllocator)
{
	_cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
}
		
void WSInversion::initializeMeasurementSet(MSProvider& msProvider, WSInversion::MSData& msData)
{
	msData.msProvider = &msProvider;
	casa::MeasurementSet& ms(msProvider.MS());
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	/**
		* Read some meta data from the measurement set
		*/
	casa::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	casa::MPosition ant1Pos = antPosColumn(0);
	
	msData.bandData = MultiBandData(ms.spectralWindow(), ms.dataDescription());
	if(Selection().HasChannelRange())
	{
		msData.startChannel = Selection().ChannelRangeStart();
		msData.endChannel = Selection().ChannelRangeEnd();
		std::cout << "Selected channels: " << msData.startChannel << '-' << msData.endChannel << '\n';
		const BandData& firstBand = msData.bandData.FirstBand();
		if(msData.startChannel >= firstBand.ChannelCount() || msData.endChannel > firstBand.ChannelCount()
			|| msData.startChannel == msData.endChannel)
		{
			std::ostringstream str;
			str << "An invalid channel range was specified! Measurement set only has " << firstBand.ChannelCount() << " channels, requested imaging range is " << msData.startChannel << " -- " << msData.endChannel << '.';
			throw std::runtime_error(str.str());
		}
	}
	else {
		msData.startChannel = 0;
		msData.endChannel = msData.bandData.FirstBand().ChannelCount();
	}
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	const MultiBandData selectedBand = msData.SelectedBand();
	if(_hasFrequencies)
	{
		_freqLow = std::min(_freqLow, selectedBand.LowestFrequency());
		_freqHigh = std::max(_freqHigh, selectedBand.HighestFrequency());
		_bandStart = std::min(_bandStart, selectedBand.BandStart());
		_bandEnd = std::max(_bandEnd, selectedBand.BandEnd());
		_startTime = std::min(_startTime, timeColumn(0).getValue().get());
	} else {
		_freqLow = selectedBand.LowestFrequency();
		_freqHigh = selectedBand.HighestFrequency();
		_bandStart = selectedBand.BandStart();
		_bandEnd = selectedBand.BandEnd();
		_startTime = timeColumn(0).getValue().get();
		_hasFrequencies = true;
	}
	
	casa::MSField fTable(ms.field());
	casa::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	casa::MDirection refDir = refDirColumn(Selection().FieldId());
	casa::MEpoch curtime = timeColumn(0);
	casa::MeasFrame frame(ant1Pos, curtime);
	casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MDirection j2000 = casa::MDirection::Convert(refDir, j2000Ref)();
	casa::Vector<casa::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];
	if(fTable.keywordSet().isDefined("WSCLEAN_DL"))
		_phaseCentreDL = fTable.keywordSet().asDouble(casa::RecordFieldId("WSCLEAN_DL"));
	else _phaseCentreDL = 0.0;
	if(fTable.keywordSet().isDefined("WSCLEAN_DM"))
		_phaseCentreDM = fTable.keywordSet().asDouble(casa::RecordFieldId("WSCLEAN_DM"));
	else _phaseCentreDM = 0.0;

	_denormalPhaseCentre = _phaseCentreDL != 0.0 || _phaseCentreDM != 0.0;
	if(_denormalPhaseCentre)
		std::cout << "Set has denormal phase centre: dl=" << _phaseCentreDL << ", dm=" << _phaseCentreDM << '\n';
	
	std::cout << "Determining min and max w & beam size... " << std::flush;
	msData.maxW = 0.0;
	msData.minW = 1e100;
	double maxBaseline = 0.0;
	std::vector<float> weightArray(selectedBand.MaxChannels());
	msProvider.Reset();
	do
	{
		size_t dataDescId;
		double uInM, vInM, wInM;
		msProvider.ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& curBand = selectedBand[dataDescId];
		double wHi = fabs(wInM / curBand.SmallestWavelength());
		double wLo = fabs(wInM / curBand.LongestWavelength());
		double baselineInM = sqrt(uInM*uInM + vInM*vInM + wInM*wInM);
		if(wHi > msData.maxW || wLo < msData.minW || baselineInM / curBand.SmallestWavelength() > maxBaseline)
		{
			msProvider.ReadWeights(weightArray.data());
			const float* weightPtr = weightArray.data();
			for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
			{
				if(*weightPtr != 0.0)
				{
					const double wavelength = curBand.ChannelWavelength(ch);
					msData.maxW = std::max(msData.maxW, fabs(wInM / wavelength));
					msData.minW = std::min(msData.minW, fabs(wInM / wavelength));
					maxBaseline = std::max(maxBaseline, baselineInM / wavelength);
				}
				++weightPtr;
			}
		}
	} while(msProvider.NextRow());
	if(msData.minW == 1e100)
	{
		msData.minW = 0.0;
		msData.maxW = 0.0;
	}
	_beamSize = 1.0 / maxBaseline;
	std::cout << "DONE (w=[" << msData.minW << " -- " << msData.maxW << "] lambdas)\n";

	if(Verbose() || !HasWGridSize())
	{
		double
			maxL = ImageWidth() * PixelSizeX() * 0.5 + fabs(_phaseCentreDL),
			maxM = ImageHeight() * PixelSizeY() * 0.5 + fabs(_phaseCentreDM),
			lmSq = maxL * maxL + maxM * maxM;
		double radiansForAllLayers;
		if(lmSq < 1.0)
			radiansForAllLayers = 2 * M_PI * (msData.maxW - msData.minW) * (1.0 - sqrt(1.0 - lmSq));
		else
			radiansForAllLayers = 2 * M_PI * (msData.maxW - msData.minW);
		size_t suggestedGridSize = size_t(ceil(radiansForAllLayers));
		if(suggestedGridSize == 0) suggestedGridSize = 1;
		if(suggestedGridSize < _cpuCount)
		{
			// When nwlayers is lower than the nr of cores, we cannot parallellize well. 
			// However, we don't want extra w-layers if we are low on mem, as that might slow down the process
			double memoryRequired = double(_cpuCount) * double(sizeof(LayeredImager::imgnum_t)) * double(ImageWidth()) * double(ImageHeight());
			if(4.0 * memoryRequired < double(_memSize))
			{
				std::cout <<
					"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
					"cores (" << _cpuCount << "). Changing suggested number of w-layers to " << _cpuCount << ".\n";
				suggestedGridSize = _cpuCount;
			}
			else {
				std::cout <<
					"The theoretically suggested number of w-layers (" << suggestedGridSize << ") is less than the number of availables\n"
					"cores (" << _cpuCount << "), but there is not enough memory available to increase the number of w-layers.\n"
					"Not all cores can be used efficiently.\n";
			}
		}
		if(Verbose())
			std::cout << "Suggested number of w-layers: " << ceil(suggestedGridSize) << '\n';
		if(!HasWGridSize())
			SetWGridSize(suggestedGridSize);
	}
}

void WSInversion::countSamplesPerLayer(MSData& msData)
{
	std::vector<size_t> sampleCount(WGridSize());
	msData.matchingRows = 0;
	msData.msProvider->Reset();
	do
	{
		double uInM, vInM, wInM;
		size_t dataDescId;
		msData.msProvider->ReadMeta(uInM, vInM, wInM, dataDescId);
		const BandData& bandData(msData.bandData[dataDescId]);
		for(size_t ch=msData.startChannel; ch!=msData.endChannel; ++ch)
		{
			double w = wInM / bandData.ChannelWavelength(ch);
			size_t wLayerIndex = _imager->WToLayer(w);
			if(wLayerIndex < WGridSize())
				++sampleCount[wLayerIndex];
		}
		++msData.matchingRows;
	} while(msData.msProvider->NextRow());
	std::cout << "Visibility count per layer: ";
	for(std::vector<size_t>::const_iterator i=sampleCount.begin(); i!=sampleCount.end(); ++i)
	{
		std::cout << *i << ' ';
	}
	std::cout << '\n';
}

void WSInversion::gridMeasurementSet(MSData &msData)
{
	const MultiBandData selectedBand(msData.SelectedBand());
	_imager->PrepareBand(selectedBand);
	std::vector<std::complex<float>> modelBuffer(selectedBand.MaxChannels());
	std::vector<float> weightBuffer(selectedBand.MaxChannels());
	
	size_t rowsRead = 0;
	msData.msProvider->Reset();
	do
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBand[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_imager->IsInLayerRange(w1, w2))
		{
			InversionWorkItem newItem;
			newItem.u = uInMeters;
			newItem.v = vInMeters;
			newItem.w = wInMeters;
			newItem.dataDescId = dataDescId;
			newItem.data = new std::complex<float>[curBand.ChannelCount()];
			
			if(DoImagePSF())
			{
				msData.msProvider->ReadWeights(newItem.data);
				if(_denormalPhaseCentre)
				{
					double shiftFactor = 2.0*M_PI* (newItem.w * (sqrt(1.0-_phaseCentreDL*_phaseCentreDL- _phaseCentreDM*_phaseCentreDM)-1.0));
					rotateVisibilities(curBand, shiftFactor, newItem.data);
				}
			}
			else {
				msData.msProvider->ReadData(newItem.data);
			}
			
			if(DoSubtractModel())
			{
				msData.msProvider->ReadModel(modelBuffer.data());
				std::complex<float>* modelIter = modelBuffer.data();
				for(std::complex<float>* iter = newItem.data; iter!=newItem.data+curBand.ChannelCount(); ++iter)
				{
					*iter -= *modelIter;
					modelIter++;
				}
			}
			msData.msProvider->ReadWeights(weightBuffer.data());
			switch(Weighting().Mode())
			{
				case WeightMode::UniformWeighted:
				case WeightMode::BriggsWeighted:
				{
					std::complex<float>* dataIter = newItem.data;
					float* weightIter = weightBuffer.data();
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						double
							u = newItem.u / curBand.ChannelWavelength(ch),
							v = newItem.v / curBand.ChannelWavelength(ch),
							weight = Weighting().IsUniform() ?
								PrecalculatedWeightInfo()->GetUniformWeight(u, v) :
								PrecalculatedWeightInfo()->GetBriggsWeight(u, v);
						*dataIter *= weight;
						_totalWeight += weight * *weightIter;
						++dataIter;
						++weightIter;
					}
				} break;
				case WeightMode::NaturalWeighted:
				{
					float* weightIter = weightBuffer.data();
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						_totalWeight += *weightIter;
						++weightIter;
					}
				} break;
				case WeightMode::DistanceWeighted:
				{
					float* weightIter = weightBuffer.data();
					double mwaWeight = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
					for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
					{
						_totalWeight += *weightIter * mwaWeight;
						++weightIter;
					}
				} break;
			}
			
			_inversionWorkLane->write(newItem);
			
			++rowsRead;
		}
	} while(msData.msProvider->NextRow());
	
	if(Verbose())
		std::cout << "Rows that were required: " << rowsRead << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsRead;
}

void WSInversion::workThreadParallel(const MultiBandData* selectedBand)
{
	std::unique_ptr<ao::lane<InversionWorkSample>[]> lanes(new ao::lane<InversionWorkSample>[_cpuCount]);
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
		lanes[i].resize(selectedBand->FirstBand().ChannelCount() * _laneBufferSize);
		bufferedLanes[i].reset(&lanes[i], std::max<size_t>(selectedBand->FirstBand().ChannelCount(), _laneBufferSize));
		
		group.add_thread(new boost::thread(&WSInversion::workThreadPerSample, this, &lanes[i]));
	}
	
	InversionWorkItem workItem;
	while(_inversionWorkLane->read(workItem))
	{
		const BandData& curBand = (*selectedBand)[workItem.dataDescId];
		InversionWorkSample sampleData;
		for(size_t ch=0; ch!=curBand.ChannelCount(); ++ch)
		{
			double wavelength = curBand.ChannelWavelength(ch);
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

void WSInversion::workThreadPerSample(ao::lane<InversionWorkSample>* workLane)
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
	msData.msProvider->ReopenRW();
	const MultiBandData selectedBandData(msData.SelectedBand());
	_imager->PrepareBand(selectedBandData);
	
	size_t rowsProcessed = 0;
	
	ao::lane<SamplingWorkItem> calcLane(_laneBufferSize+_cpuCount), writeLane(_laneBufferSize);
	boost::thread writeThread(&WSInversion::visSampleWriteThread, this, &writeLane, &msData);
	boost::thread_group calcThreads;
	for(size_t i=0; i!=_cpuCount; ++i)
		calcThreads.add_thread(new boost::thread(&WSInversion::visSampleCalcThread, this, &calcLane, &writeLane));

		
	/* Start by reading the u,v,ws in, so we don't need IO access
	 * from this thread during further processing */
	std::vector<double> us, vs, ws;
	std::vector<size_t> rowIds, dataIds;
	msData.msProvider->Reset();
	do
	{
		size_t dataDescId;
		double uInMeters, vInMeters, wInMeters;
		msData.msProvider->ReadMeta(uInMeters, vInMeters, wInMeters, dataDescId);
		const BandData& curBand(selectedBandData[dataDescId]);
		const double
			w1 = wInMeters / curBand.LongestWavelength(),
			w2 = wInMeters / curBand.SmallestWavelength();
		if(_imager->IsInLayerRange(w1, w2))
		{
			us.push_back(uInMeters);
			vs.push_back(vInMeters);
			ws.push_back(wInMeters);
			dataIds.push_back(dataDescId);
			rowIds.push_back(msData.msProvider->RowId());
			++rowsProcessed;
		}
	} while(msData.msProvider->NextRow());
	
	for(size_t i=0; i!=us.size(); ++i)
	{
		SamplingWorkItem newItem;
		newItem.u = us[i];
		newItem.v = vs[i];
		newItem.w = ws[i];
		newItem.dataDescId = dataIds[i];
		newItem.data = new std::complex<float>[selectedBandData[dataIds[i]].ChannelCount()];
		newItem.rowId = rowIds[i];
				
		calcLane.write(newItem);
	}
	std::cout << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsProcessed;
	
	calcLane.write_end();
	calcThreads.join_all();
	writeLane.write_end();
	writeThread.join();
}

void WSInversion::visSampleCalcThread(ao::lane<SamplingWorkItem>* inputLane, ao::lane<SamplingWorkItem>* outputLane)
{
	SamplingWorkItem item;
	while(inputLane->read(item))
	{
		_imager->SampleData(item.data, item.dataDescId, item.u, item.v, item.w);
		
		outputLane->write(item);
	}
}

void WSInversion::visSampleWriteThread(ao::lane<SamplingWorkItem>* samplingWorkLane, const MSData* msData)
{
	SamplingWorkItem workItem;
	while(samplingWorkLane->read(workItem))
	{
		msData->msProvider->WriteModel(workItem.rowId, workItem.data);
	}
}

void WSInversion::Invert()
{
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	_memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) _memSize / (1024.0*1024.0*1024.0);
	std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";

	MSData msDataVector[MeasurementSetCount()];
	_hasFrequencies = false;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		initializeMeasurementSet(MeasurementSet(i), msDataVector[i]);
	
	double minW = msDataVector[0].minW;
	double maxW = msDataVector[0].maxW;
	for(size_t i=1; i!=MeasurementSetCount(); ++i)
	{
		if(msDataVector[i].minW < minW) minW = msDataVector[i].minW;
		if(msDataVector[i].maxW > maxW) maxW = msDataVector[i].maxW;
	}
	
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_imager->SetGridMode(_gridMode);
	if(_denormalPhaseCentre)
		_imager->SetDenormalPhaseCentre(_phaseCentreDL, _phaseCentreDM);
	_imager->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), minW, maxW);
	_imager->SetImageImaginaryPart(ImaginaryPart());
	_imager->SetImageConjugatePart(Polarization() == Polarization::YX && ImaginaryPart());
	
	if(Verbose())
	{
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			countSamplesPerLayer(msDataVector[i]);
	}
	
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		_inversionWorkLane.reset(new ao::lane<InversionWorkItem>(16));
		
		_imager->StartInversionPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			_inversionWorkLane->clear();
			
			MSData& msData = msDataVector[i];
			
			const MultiBandData selectedBand(msData.SelectedBand());
			
			boost::thread thread(&WSInversion::workThreadParallel, this, &selectedBand);
		
			gridMeasurementSet(msData);
			
			_inversionWorkLane->write_end();
			thread.join();
		}
		
		std::cout << "Fourier transforming layers, w-term correction & summing in parallel...\n";
		_imager->FinishInversionPass();
	}
	
	if(Verbose())
	{
		size_t totalRowsRead = 0, totalMatchingRows = 0;
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
		{
			totalRowsRead += msDataVector[i].totalRowsProcessed;
			totalMatchingRows += msDataVector[i].matchingRows;
		}
		
		std::cout << "Total rows read: " << totalRowsRead << " (overhead: " << round(totalRowsRead * 100.0 / totalMatchingRows - 100.0) << "%)\n";
	}
	
	if(_denormalPhaseCentre)
		_imager->FinalizeImage(1.0/(_totalWeight * sqrt(1.0 - _phaseCentreDL*_phaseCentreDL - _phaseCentreDM*_phaseCentreDM)));
	else
		_imager->FinalizeImage(1.0/_totalWeight);
}

void WSInversion::InvertToVisibilities(const double *image)
{
	long int pageCount = sysconf(_SC_PHYS_PAGES), pageSize = sysconf(_SC_PAGE_SIZE);
	_memSize = (int64_t) pageCount * (int64_t) pageSize;
	double memSizeInGB = (double) _memSize / (1024.0*1024.0*1024.0);
	std::cout << "Detected " << round(memSizeInGB*10.0)/10.0 << " GB of system memory.\n";
	
	MSData msDataVector[MeasurementSetCount()];
	_hasFrequencies = false;
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		initializeMeasurementSet(MeasurementSet(i), msDataVector[i]);
	
	double minW = msDataVector[0].minW;
	double maxW = msDataVector[0].maxW;
	for(size_t i=1; i!=MeasurementSetCount(); ++i)
	{
		if(msDataVector[i].minW < minW) minW = msDataVector[i].minW;
		if(msDataVector[i].maxW > maxW) maxW = msDataVector[i].maxW;
	}
	
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), _cpuCount, _imageBufferAllocator, AntialiasingKernelSize(), OverSamplingFactor()));
	_imager->SetGridMode(_gridMode);
	if(_denormalPhaseCentre)
		_imager->SetDenormalPhaseCentre(_phaseCentreDL, _phaseCentreDM);
	_imager->PrepareWLayers(WGridSize(), double(_memSize)*(7.0/10.0), minW, maxW);
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		countSamplesPerLayer(msDataVector[i]);
	
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "W-term correction & fourier transforming layers to uv space (in parallel)...\n";
		if(_denormalPhaseCentre)
			_imager->PrepareImageForVisibilitySampling(image, sqrt(1.0 - _phaseCentreDL*_phaseCentreDL - _phaseCentreDM*_phaseCentreDM));
		else
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

void WSInversion::rotateVisibilities(const BandData &bandData, double shiftFactor, std::complex<float>* dataIter)
{
	for(unsigned ch=0; ch!=bandData.ChannelCount(); ++ch)
	{
		const double wShiftRad = shiftFactor / bandData.ChannelWavelength(ch);
		double rotSinD, rotCosD;
		sincos(wShiftRad, &rotSinD, &rotCosD);
		float rotSin = rotSinD, rotCos = rotCosD;
		std::complex<float> v = *dataIter;
		*dataIter = std::complex<float>(
			v.real() * rotCos  -  v.imag() * rotSin,
			v.real() * rotSin  +  v.imag() * rotCos);
		++dataIter;
	}
}
