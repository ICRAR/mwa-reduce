#include "wsinversion.h"
#include "uvwdistribution.h"
#include "imageweights.h"

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

void WSInversion::initializeMeasurementSet(const string& measurementSet, WSInversion::MSData& msData)
{
	std::cout << "Opening " << measurementSet << "... " << std::flush;
	msData.ms.reset(new casa::MeasurementSet(measurementSet));
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
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	if(_hasFrequencies)
	{
		_freqLow = std::min(_freqLow, bandData.LowestFrequency());
		_freqHigh = std::max(_freqHigh, bandData.HighestFrequency());
		_bandStart = std::min(_bandStart, bandData.BandStart());
		_bandEnd = std::max(_bandEnd, bandData.BandEnd());
		_startTime = std::min(_startTime, timeColumn(0).getValue().get());
	} else {
		_freqLow = bandData.LowestFrequency();
		_freqHigh = bandData.HighestFrequency();
		_bandStart = bandData.BandStart();
		_bandEnd = bandData.BandEnd();
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
				if(timestepIndex == IntervalStop())
				{
					msData.rowEnd = row;
					break;
				}
				time = timeColumn(row).getValue();
			}
		}
		std::cout << "DONE (" << msData.rowStart << '-' << msData.rowEnd << ")\n";
	}
	
	// Determine min and max w
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
			double wHi = fabs(wInM / bandData.SmallestWavelength());
			double wLo = fabs(wInM / bandData.LongestWavelength());
			msData.maxW = std::max(msData.maxW, wHi);
			msData.minW = std::min(msData.minW, wLo);
			maxBaseline = std::max(maxBaseline, uInM*uInM + vInM*vInM + wInM*wInM);
		}
	}
	_beamSize = bandData.SmallestWavelength() / sqrt(maxBaseline);
	std::cout << "DONE (w=[" << msData.minW << " -- " << msData.maxW << "] lambdas)\n";
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
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
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
	std::unique_ptr<casa::ROArrayColumn<casa::Complex>> modelColumn;
	if(DoSubtractModel() && !DoImagePSF())
		modelColumn.reset(new casa::ROArrayColumn<std::complex<float>>(ms, ms.columnName(casa::MSMainEnums::MODEL_DATA)));
	
	casa::IPosition dataShape = dataColumn.shape(0);
	msData.polarizationCount = dataShape[0];
	
	BandData bandData(ms.spectralWindow());
	_imager->PrepareBand(bandData);
	
	casa::Array<std::complex<float>> data(dataShape), modelData(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	size_t rowsRead = 0;
	for(size_t row=msData.rowStart; row!=msData.rowEnd; ++row)
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
				InversionWorkItem newItem;
				newItem.u = uvwArray(0);
				newItem.v = uvwArray(1);
				newItem.w = wInMeters;
				newItem.data = new std::complex<float>[msData.channelCount];
				
				weightColumn.get(row, weights);
				flagColumn.get(row, flags);
				
				double rowWeight = 1.0;
				switch(Weighting())
				{
					case UniformWeighted:
					{
						float* weightIter = weights.cbegin();
						for(size_t ch=0; ch!=bandData.ChannelCount(); ++ch)
						{
							double
								u = newItem.u / bandData.ChannelWavelength(ch),
								v = newItem.v / bandData.ChannelWavelength(ch),
								weight = PrecalculatedWeightInfo()->GetUniformWeight(u, v);
							for(size_t p=0; p!=msData.polarizationCount; ++p)
							{
								*weightIter *= weight;
								++weightIter;
							}
						}
						/*
						rowWeight = 1.0;
						double uvwDistInM = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
						float* weightIter = weights.cbegin();
						for(size_t ch=0; ch!=bandData.ChannelCount(); ++ch)
						{
							double dist = uvwDistInM / bandData.ChannelWavelength(ch); // bandData.CentreWavelength();//
							double weight = msData.uvwDistribution->WeightFromFit(dist);
							for(size_t p=0; p!=msData.polarizationCount; ++p)
							{
								*weightIter *= weight;
								++weightIter;
							}
						}*/
					} break;
					case NaturalWeighted:
						rowWeight = 1.0;
						break;
					case DistanceWeighted:
						rowWeight = sqrt(newItem.u*newItem.u + newItem.v*newItem.v + newItem.w*newItem.w);
						break;
				}
				
				dataColumn.get(row, data);
				if(DoImagePSF())
				{
					copyWeights(newItem.data, msData.channelCount, data, weights, flags, rowWeight);
				}
				else {
					if(DoSubtractModel())
					{
						modelColumn->get(row, modelData);
						casa::Array<casa::Complex>::contiter modelIter = modelData.cbegin();
						for(casa::Array<casa::Complex>::contiter iter = data.cbegin(); iter!=data.cend(); ++iter)
						{
							*iter -= *modelIter;
							modelIter++;
						}
					}
					copyWeightedData(newItem.data, msData.channelCount, data, weights, flags, rowWeight);
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
	
	BandData bandData(ms.spectralWindow());
	_imager->PrepareBand(bandData);
	
	size_t rowsProcessed = 0;
	
	_samplingWorkLane.reset(new lane<SamplingWorkItem>(16));
	boost::thread thread(&WSInversion::visSampleThread, this);
		
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
				w1 = wInMeters / bandData.LongestWavelength(),
				w2 = wInMeters / bandData.SmallestWavelength();
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
		newItem.data = new std::complex<float>[msData.channelCount];
		newItem.rowIndex = rowIndices[i];
				
		_imager->SampleData(newItem.data, newItem.u, newItem.v, newItem.w);
		
		_samplingWorkLane->write(newItem);
	}
	std::cout << "Rows that were required: " << rowsProcessed << '/' << msData.matchingRows << '\n';
	msData.totalRowsProcessed += rowsProcessed;
	
	_samplingWorkLane->write_end();
	thread.join();
		
	_modelColumn.reset();
}

void WSInversion::visSampleThread()
{
	SamplingWorkItem workItem;
	
	// Read one, which makes sure other thread has initialized the ms
	if(!_samplingWorkLane->read(workItem))
		return;
	
	casa::IPosition shape = _modelColumn->shape(0);
	casa::Array<std::complex<float>> data(shape);
	size_t channelCount = shape[1];
	
	int polIndex = polarizationIndex();
	do
	{
		_modelColumn->get(workItem.rowIndex, data);
		casa::Array<std::complex<float>>::contiter dataIter = data.cbegin();
		
		if(Polarization() == StokesI)
		{
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				if(std::isfinite(workItem.data[ch].real()))
				{
					if(AddToModel())
					{
						*dataIter += workItem.data[ch];
						*(dataIter + 3) += workItem.data[ch];
					} else {
						*dataIter = workItem.data[ch];
						*(dataIter + 3) = workItem.data[ch];
					}
				}
				dataIter += 4;
			}
		} else {
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				if(std::isfinite(workItem.data[ch].real()))
				{
					if(AddToModel())
						*(dataIter+polIndex) += workItem.data[ch];
					else
						*(dataIter+polIndex) = workItem.data[ch];
				}
				dataIter += 4;
			}
		}
		_modelColumn->put(workItem.rowIndex, data);
	}
	while(_samplingWorkLane->read(workItem));
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

	long cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), cpuCount));
	_imager->SetGridMode(_gridMode);
	_imager->PrepareWLayers(WGridSize(), memSize*2/4, minW, maxW);
	
	for(size_t i=0; i!=MeasurementSetCount(); ++i)
		countSamplesPerLayer(msDataVector[i]);
	
	_totalWeight = 0.0;
	for(size_t pass=0; pass!=_imager->NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		_inversionWorkLane.reset(new lane<InversionWorkItem>(16));
		boost::thread thread(&WSInversion::workThread, this);
		
		_imager->StartInversionPass(pass);
		
		for(size_t i=0; i!=MeasurementSetCount(); ++i)
			gridMeasurementSet(msDataVector[i]);
		
		_inversionWorkLane->write_end();
		thread.join();
		
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

	long cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	_imager = std::unique_ptr<LayeredImager>(new LayeredImager(ImageWidth(), ImageHeight(), PixelSizeX(), PixelSizeY(), cpuCount));
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

void WSInversion::copyWeightedData(std::complex<float>* dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin();
	casa::Array<float>::const_contiter weightPtr = weights.cbegin();
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin();
		
	if(Polarization() == StokesI)
	{
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			weightPtr += 3;
			inPtr += 3;
			flagPtr += 3;
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
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = *inPtr * (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			weightPtr += 4;
			inPtr += 4;
			flagPtr += 4;
		}
	}
}

void WSInversion::copyWeights(std::complex<float>* dest, size_t channelCount, const casa::Array<std::complex<float>>& data, const casa::Array<float>& weights, const casa::Array<bool>& flags, float rowWeight)
{
	casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin();
	casa::Array<float>::const_contiter weightPtr = weights.cbegin();
	casa::Array<bool>::const_contiter flagPtr = flags.cbegin();
		
	if(Polarization() == StokesI)
	{
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			}
			else {
				dest[ch] = 0;
			}
			inPtr += 3;
			weightPtr += 3;
			flagPtr += 3;
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
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			if(!*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
			{
				dest[ch] = (*weightPtr) * rowWeight;
				_totalWeight += (*weightPtr) * rowWeight;
			} else {
				dest[ch] = 0;
			}
			inPtr += 4;
			weightPtr += 4;
			flagPtr += 4;
		}
	}
}
