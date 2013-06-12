#include "banddata.h"
#include "lane.h"
#include "layeredimager.h"
#include "fitswriter.h"

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

struct WorkItem {
	double u, v, w;
	std::complex<float> *data;
};

void processWork(WorkItem &work, LayeredImager *imager)
{
	imager->AddData(work.data, work.u, work.v, work.w);
	delete[] work.data;
}

void workThread(lane<WorkItem> *workLane, LayeredImager *imager)
{
	WorkItem workItem;
	while(workLane->read(workItem))
	{
		processWork(workItem, imager);
	}
}

int main(int argc, char* argv[])
{
	int argi = 1;
	size_t imgWidth = 2048, imgHeight = 2048;
	double pixelScale = 0.01 * M_PI / 180.0;
	size_t nWLayers = 64;
	
	while(argv[argi][0] == '-')
	{
		const char *param = &argv[argi][1];
		if(strcmp(param, "size") == 0)
		{
			imgWidth = atoi(argv[argi+1]);
			imgHeight = atoi(argv[argi+2]);
			argi += 2;
		}
		else if(strcmp(param, "scale") == 0)
		{
			pixelScale = atof(argv[argi+1]) * M_PI / 180.0;
			++argi;
		}
		else if(strcmp(param, "nwlayers") == 0)
		{
			nWLayers = atoi(argv[argi+1]);
			++argi;
		}
		else {
			throw std::runtime_error("Unknown parameter");
		}
		
		++argi;
	}
	
	const char *msName(argv[argi]);
	const char *fitsfileName(argv[argi+1]);
	
	std::cout << "Opening " << msName << "... " << std::flush;
	casa::MeasurementSet ms(msName);
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
	casa::ROArrayColumn<std::complex<float> > dataColumn(ms, "DATA");
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
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
	double phaseCentreRA = j2000Val[0];
	double phaseCentreDec = j2000Val[1];
	
	std::cout << 'D' << std::flush;
	casa::IPosition dataShape = dataColumn.shape(0);
	casa::Array<std::complex<float> > data(dataShape);
	casa::Array<bool> flags(dataShape);
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

	LayeredImager imager(imgWidth, imgHeight, pixelScale);
	imager.PrepareForObservation(nWLayers, memSize*2/4, minW, maxW, bandData);
	
	std::vector<size_t> sampleCount(nWLayers);
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
				++sampleCount[imager.WToLayer(w)];
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
	
	size_t totalRowsRead = 0;
	for(size_t pass=0; pass!=imager.NPasses(); ++pass)
	{
		std::cout << "Starting gridding pass " << pass << ".\n";
		lane<WorkItem> workLane(16);
		boost::thread thread(&workThread, &workLane, &imager);
		
		imager.StartPass(pass);
		
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
				if(imager.IsInLayerRange(w1, w2))
				{
					dataColumn.get(row, data);
					flagColumn.get(row, flags);
							
					WorkItem newItem;
					newItem.data = new std::complex<float>[channelCount];
					newItem.u = uvwArray(0);
					newItem.v = uvwArray(1);
					newItem.w = wInMeters;
					casa::Array<std::complex<float> >::const_contiter inPtr = data.cbegin();
					casa::Array<bool>::const_contiter flagPtr = flags.cbegin();
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						if(*flagPtr)
							newItem.data[ch] = 0;
						else
							newItem.data[ch] = *inPtr; // copy XX for now
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							++inPtr;
							++flagPtr;
						}
					}
					workLane.write(newItem);
					
					++rowsRead;
				}
			}
		}
		std::cout << "Pass " << pass << ", rows that were required: " << rowsRead << '/' << matchingRows << '\n';
		totalRowsRead += rowsRead;
		
		workLane.write_end();
		thread.join();
		
		std::cout << "Summing down layers...\n";
		imager.FinishPass();
	}
	std::cout << "Total rows read: " << totalRowsRead << " (overhead: " << round(totalRowsRead * 100.0 / matchingRows - 100.0) << "%)\n";
	
	std::cout << "Writing image... " << std::flush;
	imager.FinalizeImage();
	FitsWriter writer(fitsfileName);
	writer.Write(imager.Image(), imgWidth, imgHeight, phaseCentreRA, phaseCentreDec, -pixelScale, pixelScale);
	std::cout << "DONE\n";
}
