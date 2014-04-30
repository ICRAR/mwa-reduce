#include "beam/tilebeam.h"

#include "fitsreader.h"
#include "imagecoordinates.h"
#include "fitswriter.h"
#include "banddata.h"
#include "matrix2x2.h"
#include "numberlist.h"
#include "progressbar.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>

#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <stdexcept>

template<typename TileImplementation, bool DoSquare>
void MakeBeam(double** imgPtr, size_t width, size_t height, double pixelSizeX, double pixelSizeY, double refRA, double refDec, double arrLatitude, double zenithHa, double zenithDec, double centralFrequency, const double* delays, const casa::MeasFrame& frame, double dl, double dm)
{
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MDirection::Convert
		j2000ToHaDecRef(j2000Ref, hadecRef),
		j2000ToAzelGeoRef(j2000Ref, azelgeoRef);

	TileBeamBase<TileImplementation> tilebeam(delays);
	ProgressBar progressBar("Constructing beam");
	for(size_t y=0;y!=height;++y)
	{
		for(size_t x=0;x!=width;++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, pixelSizeX, pixelSizeY, width, height, l, m);
			l += dl; m += dm;
			ImageCoordinates::LMToRaDec(l, m, refRA, refDec, ra, dec);
			
			std::complex<double> gain[4];
			tilebeam.ArrayResponse(ra, dec, j2000Ref, j2000ToHaDecRef, j2000ToAzelGeoRef, arrLatitude, zenithHa, zenithDec, centralFrequency, gain);
			if(DoSquare) {
				std::complex<double> gainSq[4];
				Matrix2x2::ATimesHermB(gainSq, gain, gain);
				Matrix2x2::Assign(gain, gainSq);
			}
			
			for(size_t i=0; i!=4; ++i)
			{
				*imgPtr[i*2] = gain[i].real();
				*imgPtr[i*2 + 1] = gain[i].imag();
				++imgPtr[i*2];
				++imgPtr[i*2 + 1];
			}
		}
		progressBar.SetProgress(y, height);
	}
	progressBar.SetProgress(height, height);
}

int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		std::cout << "Syntax: beam [-2013 | -2014] [-allsky] [-square] [-proto <input fitsfile>] [-name <prefix>] [-ms <measurementset>] [-delays <0,0,..>]\n";
		return 0;
	}
	
	int argi= 1;
	const char* inpFitsname = 0;
	const char* msName = 0;
	string prefixName = "beam"; 
	double delays[16];
	bool doSquare = false;
	bool use2013 = true;
	for(size_t i=0; i!=16; ++i)
		delays[i] = 0.0;
	while(argi<argc && argv[argi][0] == '-')
	{
		std::string param(&argv[argi][1]);
		if(param == "proto")
		{
			++argi;
			inpFitsname = argv[argi];
		}
		else if(param == "ms")
		{
			++argi;
			msName = argv[argi];
		}
		else if(param == "name")
		{
			++argi;
			prefixName = argv[argi];
		}
		else if(param == "2013")
			use2013 = true;
		else if(param == "2014")
			use2013 = false;
		else if(param == "allsky")
		{
		}
		else if(param == "delays")
		{
			++argi;
			std::vector<int> list;
			NumberList::ParseIntList(argv[argi], list);
			if(list.size() != 16) {
				std::cerr << "Need 16 delays\n";
			  exit(1);
			}
			for(size_t i=0; i!=16; ++i)
				delays[i] = list[i];
		}
		else if(param == "square")
		{
			doSquare = true;
		}
		else throw std::runtime_error(std::string("Invalid param: ") + param);
		++argi;
	}
	
	casa::MPosition arrayPos;
	casa::MEpoch time;
	double centralFrequency;
	
	if(msName == 0)
	{
		arrayPos = casa::MPosition(casa::MVPosition(-2.55952e+06, 5.09585e+06, -2.84899e+06)); // pos of tile 011
		time = casa::MEpoch(casa::MVEpoch(casa::Quantity(4.88193e+09, "s")));
		centralFrequency = 150000000.0;
	}
	else {
		/**
			* Read some meta data from the measurement set
			*/
		casa::MeasurementSet ms(msName);
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::MSAntenna aTable = ms.antenna();
		if(aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
		casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
		arrayPos = antPosColumn(0);
		casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		// Find the mid time step
		casa::MEpoch firstTime = timeColumn(0);
		casa::MEpoch lastTime = timeColumn(ms.nrow()-1);
		std::cout << "Start time = " << firstTime << ", mid time = " << lastTime << '\n';
		time = casa::MEpoch(casa::MVEpoch(0.5 * (firstTime.getValue().get() + lastTime.getValue().get())), firstTime.getRef());
		std::cout << "Using mid time = " << time << '\n';
		
		casa::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
		casa::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
		casa::Array<int> delaysArr = delaysCol(0);
		casa::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
		for(int i=0; i!=16; ++i)
			delays[i] = delaysArrPtr[i];
		
		BandData bandData(ms.spectralWindow());
		centralFrequency = bandData.CentreFrequency();
	}
	
	std::cout << "Delays: [";
	for(int i=0; i!=16; ++i)
	{
		std::cout << delays[i];
		if(i != 15) std::cout << ',';
	}
	std::cout << "]\n";
	
	casa::MeasFrame frame(arrayPos, time);
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MPosition wgs = casa::MPosition::Convert(arrayPos, casa::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // arrayPos.getValue().getLat();
	
	casa::MDirection zenith(casa::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
	casa::MDirection zenithHaDec = casa::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getAngle().getValue()[0];
	double zenithDec = zenithHaDec.getAngle().getValue()[1];
	std::cout << "Zenith: "
		<< (casa::MDirection::Convert(zenith, j2000Ref)()).getAngle().getValue()[0]*180.0/M_PI << " RA, "
		<< zenithDec*180.0/M_PI << " dec, "
		<< zenithHa*180.0/M_PI << " HA.\n";
		
	size_t width, height;
	double pixelSizeX, pixelSizeY, refRA, refDec, phaseCentreDL, phaseCentreDM;
	FitsWriter writer;
	if(inpFitsname == 0)
	{
		// All sky
		width = 512;
		height = 512;
		pixelSizeX = 2.0 / (double) width;
		pixelSizeY = 2.0 / (double) height;
		refRA = (casa::MDirection::Convert(zenith, j2000Ref)()).getAngle().getValue()[0];
		refDec = zenithDec;
		phaseCentreDL = 0.0;
		phaseCentreDM = 0.0;
		double bandWidth = 1000000.0;
		
		writer.SetImageDimensions(width, height, refRA, refDec, pixelSizeX, pixelSizeY);
		writer.SetFrequency(centralFrequency, bandWidth);
	}
	else {
		FitsReader reader(inpFitsname);
		writer.SetMetadata(reader);
		
		width = reader.ImageWidth();
		height = reader.ImageHeight();
		pixelSizeX = reader.PixelSizeX();
		pixelSizeY = reader.PixelSizeY();
		refRA = reader.PhaseCentreRA();
		refDec = reader.PhaseCentreDec();
		phaseCentreDL = reader.PhaseCentreDL();
		phaseCentreDM = reader.PhaseCentreDM();
	}
	std::cout << "Reference dir: "
		<< refRA*180.0/M_PI << " RA, "
		<< refDec*180.0/M_PI << " dec.\n";
	
	
	std::vector<double> outImage[8];
	double *imgPtr[8];
	for(size_t i=0; i!=8; ++i)
	{
		outImage[i].resize(width*height);
		imgPtr[i] = &outImage[i][0];
	}

	if(doSquare)
	{
		if(use2013)
			MakeBeam<TileBeam2013,true>(imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec, arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame, phaseCentreDL, phaseCentreDM);
		else
			MakeBeam<TileBeam2014,true>(imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec, arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame, phaseCentreDL, phaseCentreDM);
	}
	else {
		if(use2013)
			MakeBeam<TileBeam2013,false>(imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec, arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame, phaseCentreDL, phaseCentreDM);
		else
			MakeBeam<TileBeam2014,false>(imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec, arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame, phaseCentreDL, phaseCentreDM);
	}
	
	std::cout << "\nWriting...\n";
	
	const string names[8] = {
		prefixName+"-xxr.fits", prefixName+"-xxi.fits", prefixName+"-xyr.fits", prefixName+"-xyi.fits", 
		prefixName+"-yxr.fits", prefixName+"-yxi.fits", prefixName+"-yyr.fits", prefixName+"-yyi.fits"
	};
	
	for(size_t i=0; i!=8; ++i)
		writer.Write<double>(names[i], &outImage[i][0]);
}
