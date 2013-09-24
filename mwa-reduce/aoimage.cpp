#include "btpimager.h"
#include "lane.h"
#include "banddata.h"
#include "fitswriter.h"
#include "imageweights.h"
#include "model.h"
#include "modelrenderer.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MCDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <stdexcept>

#include <unistd.h>

#include <boost/thread/thread.hpp>

using namespace casa;

typedef BTPImager::NumType NumType;
typedef BTPImager::ImageNum ImageNum;

struct ImageWork
{
	NumType uTimesLambda, vTimesLambda, wTimesLambda;
	NumType zenithDistance, paralacticAngle, weight;
	std::complex<float> *data;
};

struct ThreadParameters
{
	lane<ImageWork> *worklane;
	BTPImager *imager;
	size_t threadIndex;
};

void workFunction(const ThreadParameters parameters)
{
	lane<ImageWork> &worklane = *parameters.worklane;
	BTPImager &imager = *parameters.imager;
	size_t threadIndex = parameters.threadIndex;
	ImageWork work;
	while(worklane.read(work))
	{
		if(work.weight != 0.0)
			imager.AddTimestep(threadIndex, work.uTimesLambda, work.vTimesLambda, work.wTimesLambda, work.data, work.zenithDistance, work.paralacticAngle, work.weight);
		delete[] work.data;
	}
}

enum Polarization { StokesIPol, XXPol, YYPol, PsfPol };

struct ImageInfo
{
	double phaseCentreRA, phaseCentreDec;
	double highestFrequency, lowestFrequency, bandwidth, dateObs;
	bool onlyModel, haveTimeRange, haveUVRange;
	size_t timeRangeStart, timeRangeStop;
	size_t uvRangeStart, uvRangeStop;
	enum Polarization polarization;
};

size_t readDataStokesI(size_t channelCount, size_t polarizationCount,
								std::complex<float> *outPtr,
								bool *outFlagPtr,
								casa::Array<std::complex<float> >::const_iterator inPtr,
								casa::Array<bool>::const_iterator flagPtr)
{
	size_t sampleCount = 0;
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		bool hasSample = false;
		if(polarizationCount == 1)
		{
			if(*flagPtr)
				*outPtr = std::complex<float>(0.0, 0.0);
			else {
				*outPtr = std::complex<float>(inPtr->real(), inPtr->imag());
				hasSample = true;
			}
			++inPtr;
			++flagPtr;
		} else if(polarizationCount == 2)
		{
			bool flagXX = *flagPtr;
			++flagPtr;
			bool flagYY = *flagPtr;
			++flagPtr;
			
			if(flagXX || flagYY) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr;
				++inPtr;
			}
			else {
				std::complex<float> xx = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;
				*outPtr = std::complex<float>((xx.real() + inPtr->real())*0.5, (xx.imag() + inPtr->imag())*0.5);
				++inPtr;
				hasSample = true;
			}
		} else if(polarizationCount == 4)
		{
			bool flagXX = *flagPtr;
			++flagPtr;++flagPtr;++flagPtr;
			bool flagYY = *flagPtr;
			++flagPtr;
			if(flagXX || flagYY) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr; ++inPtr; ++inPtr; ++inPtr;
			}
			else {
				std::complex<float> xx = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;++inPtr;++inPtr;
				*outPtr = std::complex<float>((xx.real() + inPtr->real())*0.5, (xx.imag() + inPtr->imag())*0.5);
				++inPtr;
				hasSample = true;
			}
		}
		if(!std::isfinite(outPtr->real()) || !std::isfinite(outPtr->imag()))
		{
			*outPtr = std::complex<float>(0.0, 0.0);
			hasSample = false;
		}
		else if(hasSample)
			sampleCount++;
		*outFlagPtr = !hasSample;
		++outFlagPtr;
		++outPtr;
	}
	return sampleCount;
}

size_t readDataXX(size_t channelCount, size_t polarizationCount,
								std::complex<float> *outPtr,
								bool *outFlagPtr,
								casa::Array<std::complex<float> >::const_iterator inPtr,
								casa::Array<bool>::const_iterator flagPtr)
{
	size_t sampleCount = 0;
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		bool hasSample = false;
		if(polarizationCount == 2)
		{
			bool flagXX = *flagPtr;
			++flagPtr;
			++flagPtr;
			
			if(flagXX) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr;
				++inPtr;
			}
			else {
				*outPtr = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;
				++inPtr;
				hasSample = true;
			}
		} else if(polarizationCount == 4)
		{
			bool flagXX = *flagPtr;
			++flagPtr;++flagPtr;++flagPtr;++flagPtr;
			if(flagXX) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr; ++inPtr; ++inPtr; ++inPtr;
			}
			else {
				*outPtr = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;++inPtr;++inPtr;++inPtr;
				hasSample = true;
			}
		}
		if(!std::isfinite(outPtr->real()) || !std::isfinite(outPtr->imag()))
		{
			*outPtr = std::complex<float>(0.0, 0.0);
			hasSample = false;
		}
		else if(hasSample)
			sampleCount++;
		*outFlagPtr = !hasSample;
		++outFlagPtr;
		++outPtr;
	}
	return sampleCount;
}

size_t readDataYY(size_t channelCount, size_t polarizationCount,
								std::complex<float> *outPtr,
								bool *outFlagPtr,
								casa::Array<std::complex<float> >::const_iterator inPtr,
								casa::Array<bool>::const_iterator flagPtr)
{
	size_t sampleCount = 0;
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		bool hasSample = false;
		if(polarizationCount == 2)
		{
			++flagPtr;
			bool flagYY = *flagPtr;
			++flagPtr;
			
			if(flagYY) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr;
				++inPtr;
			}
			else {
				++inPtr;
				*outPtr = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;
				hasSample = true;
			}
		} else if(polarizationCount == 4)
		{
			++flagPtr;++flagPtr;++flagPtr;
			bool flagYY = *flagPtr;
			++flagPtr;
			if(flagYY) {
				*outPtr = std::complex<float>(0.0, 0.0);
				++inPtr; ++inPtr; ++inPtr; ++inPtr;
			}
			else {
				++inPtr;++inPtr;++inPtr;
				*outPtr = std::complex<float>(inPtr->real(), inPtr->imag());
				++inPtr;
				hasSample = true;
			}
		}
		if(!std::isfinite(outPtr->real()) || !std::isfinite(outPtr->imag()))
		{
			*outPtr = std::complex<float>(0.0, 0.0);
			hasSample = false;
		}
		else if(hasSample)
			sampleCount++;
		*outFlagPtr = !hasSample;
		++outFlagPtr;
		++outPtr;
	}
	return sampleCount;
}

size_t readDataWeights(size_t channelCount, size_t polarizationCount,
								std::complex<float> *outPtr,
								bool *outFlagPtr,
								casa::Array<std::complex<float> >::const_iterator inPtr,
								casa::Array<bool>::const_iterator flagPtr)
{
	size_t sampleCount = 0;
	for(size_t ch=0;ch!=channelCount;++ch)
	{
		bool hasSample = false;
		if(polarizationCount == 1)
		{
			if(*flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag()))
				hasSample = true;
			++inPtr;
			++flagPtr;
		} else if(polarizationCount == 2)
		{
			bool flagXX = *flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag());
			++flagPtr; ++inPtr;
			bool flagYY = *flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag());
			++flagPtr; ++inPtr;
			
			hasSample = !(flagXX || flagYY);
		} else if(polarizationCount == 4)
		{
			bool flagXX = *flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag());
			++flagPtr;++flagPtr;++flagPtr;
			++inPtr; ++inPtr; ++inPtr;
			bool flagYY = *flagPtr && std::isfinite(inPtr->real()) && std::isfinite(inPtr->imag());
			++flagPtr;
			++inPtr;
			hasSample = !(flagXX || flagYY);
		}
		if(hasSample)
			++sampleCount;
		*outFlagPtr = !hasSample;
		*outPtr = hasSample ? 1.0 : 0.0;
		++outFlagPtr;
		++outPtr;
	}
	return sampleCount;
}

size_t readData(enum Polarization polarization,
								size_t channelCount, size_t polarizationCount,
								std::complex<float> *outPtr,
								bool *outFlagPtr,
								casa::Array<std::complex<float> >::const_iterator inPtr,
								casa::Array<bool>::const_iterator flagPtr)
{
	switch(polarization)
	{
		case StokesIPol:
			return readDataStokesI(channelCount, polarizationCount,
				outPtr, outFlagPtr, inPtr, flagPtr);
		case XXPol:
			return readDataXX(channelCount, polarizationCount,
				outPtr, outFlagPtr, inPtr, flagPtr);
		case YYPol:
			return readDataYY(channelCount, polarizationCount,
				outPtr, outFlagPtr, inPtr, flagPtr);
		case PsfPol:
			return readDataWeights(channelCount, polarizationCount,
				outPtr, outFlagPtr, inPtr, flagPtr);
	}
	throw std::runtime_error("Unsupported polarization");
}

void image(const char *msName, const char *columnName, BTPImager &imager, size_t avgFactor, ImageInfo &info)
{
	std::cout << "Opening " << msName << "... " << std::flush;
	MeasurementSet ms(msName);
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	/**
		* Read some meta data from the measurement set
		*/
	std::cout << 'A' << std::flush;
	MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(MSAntennaEnums::POSITION));
	MPosition ant1Pos = antPosColumn(0);
	
	std::cout << 'B' << std::flush;
	BandData bandData(ms.spectralWindow());
	size_t channelCount = bandData.ChannelCount();
	
	std::cout << 'F' << std::flush;
	MSField fTable(ms.field());
	if(fTable.nrow() != 1) throw std::runtime_error("Need exactly one field in set");
	MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(MSFieldEnums::REFERENCE_DIR));
	MDirection refDir = refDirColumn(0);
	
	std::cout << 'C' << std::flush;
	typedef float num_t;
	typedef std::complex<num_t> complex_t;
	ROScalarColumn<int> ant1Column(ms, ms.columnName(MSMainEnums::ANTENNA1));
	ROScalarColumn<int> ant2Column(ms, ms.columnName(MSMainEnums::ANTENNA2));
	ROArrayColumn<std::complex<float> > dataColumn(ms, columnName);
	ROArrayColumn<bool> flagColumn(ms, ms.columnName(MSMainEnums::FLAG));
	ROArrayColumn<double> uvwColumn(ms, ms.columnName(MSMainEnums::UVW));
	MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(MSMainEnums::TIME));
	
	std::cout << 'D' << std::flush;
	IPosition dataShape = dataColumn.shape(0);
	Array<std::complex<float> > data(dataShape);
	Array<bool> flags(dataShape);
	unsigned polarizationCount = dataShape[0];
	std::cout << " DONE\n";
	
	NumType
		highestFrequency = (bandData.HighestFrequency()*2.0 - bandData.FrequencyStep()*(avgFactor-1))*0.5,
		frequencyStep = bandData.FrequencyStep()*avgFactor;
	std::cout << "ChannelCount: " << channelCount << ", polarizationCount: " << polarizationCount << ", freqstep: " << frequencyStep << '\n';
	
	std::cout << "Initializing weights... " << std::flush;
	ImageWeights weights(imager.ImageSize(), imager.ImageSize(), imager.PixelScale());
	bool weightsNeedData = true;
	if(info.onlyModel) weightsNeedData = false;
	
	double minDist = 1e100, maxDist = 0.0;
	for(size_t row=0;row!=ms.nrow();++row)
	{
		if(ant1Column(row) != ant2Column(row))
		{
			std::complex<float> formattedData[channelCount];
			bool formattedFlags[channelCount];
			casa::Array<double> uvwArray = uvwColumn(row);
			casa::Array<double>::const_iterator iter = uvwArray.begin();
			double u = *iter;
			++iter;
			double v = *iter;
			double dist = sqrt(u*u + v*v);
			if(dist < minDist) minDist = dist;
			if(dist > maxDist) maxDist = dist;
			
			if(weightsNeedData)
			{
				dataColumn.get(row, data);
				flagColumn.get(row, flags);
				readData(info.polarization, channelCount, polarizationCount, formattedData, formattedFlags, data.begin(), flags.begin());
				weights.Grid(formattedData, formattedFlags, u, v, channelCount, highestFrequency-frequencyStep*channelCount/avgFactor, frequencyStep);
			}
		}
	}
	std::cout << "DONE (min,max uv=" << minDist << ',' << maxDist << "m)\n";

	std::cout << "Initializing imager...\n";
	imager.Initialize(minDist, maxDist, highestFrequency, frequencyStep, channelCount/avgFactor);
	std::cout << "Done initializing imager.\n";
	
	MPosition wgs = MPosition::Convert(ant1Pos, MPosition::WGS84)();
	//use wgs or ant1Pos?! differ by 1%
	double latitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	MEpoch curtime = timeColumn(0);
	size_t curTimeIndex = (size_t) (-1);
	MeasFrame frame(ant1Pos, curtime);
	MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	MDirection hadec = MDirection::Convert(refDir, hadecRef)();
	MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	MDirection j2000 = MDirection::Convert(refDir, j2000Ref)();
	curtime = MEpoch(curtime.getValue() - 1, curtime.getRef()); // trigger calculation of pAngle, etc.
		
	if(!info.onlyModel)
	{
		/**
		* Start threads
		*/
		size_t cpuCount = imager.ImageCount();
		std::cout << "Starting " << cpuCount << " work threads... " << std::flush;
		
		lane<ImageWork> worklane(cpuCount * 2);
		ThreadParameters parameters;
		parameters.worklane = &worklane;
		parameters.imager = &imager;
		boost::thread *threads[cpuCount];
		for(size_t i=0;i!=cpuCount;++i)
		{
			parameters.threadIndex = i;
			threads[i] = new boost::thread(&workFunction, parameters);
		}
		std::cout << "DONE\n";
		
		/**
		* Push all work
		*/
		double zenithDistance, pAngle;
		bool formattedFlags[channelCount];
		std::cout << "Making image";
		for(size_t row=0;row!=ms.nrow();++row)
		{
			if(ant1Column(row) != ant2Column(row))
			{
				bool isSelected = true;
				
				if(curtime.getValue() != timeColumn(row).getValue())
				{
					curtime = timeColumn(row);
					++curTimeIndex;
					frame.set(curtime);
					MDirection::Ref ref(casa::MDirection::HADEC, frame);
					hadec = MDirection::Convert(refDir, ref)();
					//std::cout << refDir << ", HA, DEC of phasedir=" << hadec << ',' << hadec.getValue() << '\n';
					Vector<Double> hadecVal = hadec.getValue().get();
					double ha = hadecVal[0];
					double dec = hadecVal[1];
					double sinLat = sin(latitude), cosLat = cos(latitude);
					double sinDec = sin(dec), cosDec = cos(dec);
					double cosHA = cos(ha), sinHA = sin(ha);
					zenithDistance = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
					//This is Perley's version:
					pAngle = atan2(cosLat * sinHA, sinLat * cosDec - cosLat * sinDec * cosHA);
					//std::cout << "HA=" << (ha*180/M_PI) << ", DEC=" << (dec*180/M_PI) << ", zenith dist=acos(" << (sinLat * sinDec + cosLat * cosDec * cosHA) << ")=" << (zenithDistance*180/M_PI) << ", paralactic angle=" << (pAngle*180/M_PI) << '\n';
					
					MDirection::Ref ref2(casa::MDirection::AZELGEO, frame);
					MDirection azel = MDirection::Convert(refDir, ref2)();
					Vector<Double> azelVal = azel.getValue().get();
					double azimuth = azelVal[0];
					double zenithDistance2 = asin(cosDec * sinHA / sin(azimuth));
					double pAngle3 = asin(-sin(azimuth) * cosLat / cosDec);
					//std::cout << "Z=" << (zenithDistance*180/M_PI) << ", Z2=" << (zenithDistance2*180/M_PI) << ", pangle=" << (pAngle*180.0/M_PI) << ", pangleold=" << (atan2(cosLat * sinHA, sinLat * cosDec - cosLat * sinDec * cosHA)*180/M_PI) << "pAngle3=" << (pAngle3*180.0/M_PI) << '\n';
					
					pAngle = -pAngle;
					// I think this is because of southern hemisphere.
					// For Northern hemisphere, no negation is needed.
					zenithDistance = -zenithDistance;
					
					std::cout << '.' << std::flush;
				}
				
				if(info.haveTimeRange)
				{
					if(curTimeIndex < info.timeRangeStart || curTimeIndex >= info.timeRangeStop)
						isSelected = false;
				}
				
				ImageWork work;
				casa::Vector<double> uvw = uvwColumn(row);
				casa::Vector<double>::const_iterator uvwIter = uvw.begin();
				work.uTimesLambda = *uvwIter; ++uvwIter;
				work.vTimesLambda = *uvwIter; ++uvwIter;
				work.wTimesLambda = *uvwIter;
				
				if(info.haveUVRange)
				{
					double uv = sqrt(work.uTimesLambda*work.uTimesLambda + 
						work.vTimesLambda*work.vTimesLambda);
					if(uv < info.uvRangeStart || uv > info.uvRangeStop)
						isSelected = false;
				}
				
				if(isSelected)
				{
					dataColumn.get(row, data);
					flagColumn.get(row, flags);
					
					work.zenithDistance = zenithDistance;
					work.paralacticAngle = pAngle;
					std::complex<float> *outPtr = new std::complex<float>[channelCount];
					work.data = outPtr;
					casa::Array<std::complex<float> >::const_iterator inPtr = data.begin();
					casa::Array<bool>::const_iterator flagPtr = flags.begin();
					
					size_t sampleCount = readData(info.polarization, channelCount, polarizationCount, outPtr, formattedFlags, inPtr, flagPtr);
					
					if(avgFactor != 1)
					{
						outPtr = work.data;
						std::complex<float> *avgInPtr = work.data;
						for(size_t ch=0;ch!=channelCount/avgFactor;++ch)
						{
							std::complex<float> avgSample(0.0, 0.0);
							for(size_t a=0;a!=avgFactor;++a)
							{
								avgSample += *avgInPtr;
								++avgInPtr;
							}
							*outPtr = std::complex<float>(avgSample.real() / avgFactor, avgSample.imag() / avgFactor);
							++outPtr;
						}
					}
					double weight = weights.ApplyWeights(work.data, formattedFlags, work.uTimesLambda, work.vTimesLambda, channelCount, highestFrequency-frequencyStep*channelCount/avgFactor, frequencyStep);
					work.weight = weight * avgFactor;
					worklane.write(work);
				}
			}
		}
		std::cout << '\n';
		worklane.write_end();
		
		for(size_t i=0;i!=cpuCount;++i)
		{
			std::cout << "Waiting for thread " << i << " to finish... " << std::flush;
			threads[i]->join();
			delete threads[i];
			std::cout << "DONE\n";
		}
	}
	Vector<Double> j2000Val = j2000.getValue().get();
	info.phaseCentreRA = j2000Val[0];
	info.phaseCentreDec = j2000Val[1];
	info.highestFrequency = std::max(info.highestFrequency, bandData.HighestFrequency());
	info.lowestFrequency = std::min(info.lowestFrequency, bandData.LowestFrequency());
	info.bandwidth = std::max(bandData.Bandwidth(), info.highestFrequency - info.lowestFrequency);
	info.dateObs = timeColumn(0).getValue().get();
}

int main(int argc, char *argv[])
{
	std::cout <<
		"Welcome to the AOImager. This imager will perform an exact inversion of the\n"
		"visibility function. It applies exact widefield corrections, uses frequency\n"
		"low-pass filtering techniques to reduce sidelobes from off-axis sources and\n"
		"applies quasi non-periodic boundary conditions to avoid aliasing.\n\n";
	
	int argi = 1;
	
	size_t avgFactor = 1;
	NumType pixelScale = 0.1*(M_PI/180.0); // '0.05 deg'
	const char *columnName = "DATA", *modelFilename = 0;
	bool onlyModel = false;
	enum Polarization pol = StokesIPol;
	bool haveTimeRange = false, haveUVRange = false;
	size_t timeRangeStart, timeRangeStop;
	size_t uvRangeStart, uvRangeStop;
	while(argc - argi >= 1 && argv[argi][0]=='-')
	{
		if(argc - argi >= 2 && strcmp(argv[argi], "-avg")==0)
		{
			++argi;
			avgFactor = atoi(argv[argi]);
		}
		else if(argc - argi >= 2 && strcmp(argv[argi], "-scale")==0)
		{
			++argi;
			pixelScale = atof(argv[argi])*(M_PI/180.0);
		}
		else if(argc - argi >= 2 && strcmp(argv[argi], "-column")==0)
		{
			++argi;
			columnName = argv[argi];
		}
		else if(argc - argi >= 2 && strcmp(argv[argi], "-m")==0)
		{
			onlyModel = true;
		}
		else if(argc - argi >= 2 && strcmp(argv[argi], "-model")==0)
		{
			++argi;
			modelFilename = argv[argi];
		}
		else if(argc - argi >= 3 && strcmp(argv[argi], "-timerange")==0)
		{
			haveTimeRange = true;
			++argi;
			timeRangeStart = atoi(argv[argi]);
			++argi;
			timeRangeStop = atoi(argv[argi]);
		}
		else if(argc - argi >= 3 && strcmp(argv[argi], "-uvrange")==0)
		{
			haveUVRange = true;
			++argi;
			uvRangeStart = atoi(argv[argi]);
			++argi;
			uvRangeStop = atoi(argv[argi]);
		}
		else if(strcmp(argv[argi], "-xx")==0)
		{
			pol = XXPol;
		}
		else if(strcmp(argv[argi], "-yy")==0)
		{
			pol = YYPol;
		}
		else if(strcmp(argv[argi], "-psf")==0)
		{
			pol = PsfPol;
		}
		else throw std::runtime_error(std::string("Unknown parameter or incorrectly used: ") + argv[argi]);
		
		++argi;
	}
	if(argc - argi < 3)
	{
		std::cout << "Syntax: aoimager [-avg <num>] [-scale <pixelscale in deg>] [-column <name>] [-m] [-model <file>] <imgsize> <fitsfile> <measurement set> [..more measurement sets..]\n";
		return -1;
	}
	const size_t imgSize = atoi(argv[argi]);
	const char *fitsfile = argv[argi+1];
	argi+=2;
	
	std::cout << "Please stand by while I'll prepare " << fitsfile << " for you...\n\n";
	
	std::cout << "Constructing imager... " << std::flush;
	long cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	BTPImager imager(cpuCount, imgSize, pixelScale);
	ImageInfo imageInfo;
	imageInfo.lowestFrequency = 1e30;
	imageInfo.highestFrequency = 0.0;
	imageInfo.bandwidth = 0.0;
	imageInfo.dateObs = 0.0;
	imageInfo.onlyModel = onlyModel;
	imageInfo.polarization = pol;
	imageInfo.haveTimeRange = haveTimeRange;
	imageInfo.timeRangeStart = timeRangeStart;
	imageInfo.timeRangeStop = timeRangeStop;
	imageInfo.haveUVRange = haveUVRange;
	imageInfo.uvRangeStart = uvRangeStart;
	imageInfo.uvRangeStop = uvRangeStop;
	std::cout << "DONE\n";
	
	while(argi < argc)
	{
		const char *msName = argv[argi];
	
		image(msName, columnName, imager, avgFactor, imageInfo);
		
		ImageNum *imageData = new ImageNum[imager.ImageSize() * imager.ImageSize()];
		imager.GetIntermediateResult(imageData);
		
		if(imager.SkippedTimesteps() != 0)
			std::cout << "Skipped " << imager.SkippedTimesteps() << " timesteps because their frequency was higher than imaging resolution.\n";
		
		if(modelFilename != 0)
		{
			Model model(modelFilename);
			std::cout << "Rendering " << model.SourceCount() << " sources to image... " << std::flush;
			ModelRenderer renderer(imageInfo.phaseCentreRA, imageInfo.phaseCentreDec, pixelScale, pixelScale);
			renderer.Restore(imageData, imager.ImageSize(), imager.ImageSize(), model, 0.5 / imager.OverallMaxUVDist(), imageInfo.lowestFrequency, imageInfo.highestFrequency, 0);
			std::cout << "DONE\n";
		}
		
		++argi;
		if(argi < argc)
			std::cout << "Writing intermediate fits file... " << std::flush;
		else
			std::cout << "Writing final fits file... " << std::flush;
		FitsWriter writer(fitsfile);
		writer.Write(imageData, imager.ImageSize(), imager.ImageSize(), imageInfo.phaseCentreRA, imageInfo.phaseCentreDec, pixelScale, pixelScale, (imageInfo.lowestFrequency + imageInfo.highestFrequency) * 0.5, imageInfo.bandwidth, imageInfo.dateObs);
		delete[] imageData;
		std::cout << "DONE\n";
	}
}
