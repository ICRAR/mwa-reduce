#include "tilebeam.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "fitswriter.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MCPosition.h>

#include <measures/TableMeasures/ScalarMeasColumn.h>

#include <stdexcept>

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cout << "Syntax: beam <input fitsfile> <frequency-in-hz> <measurementset>\n";
		return 0;
	}
	const char *inpFitsname = argv[1];
	double frequency = atoi(argv[2]);
	const char *msName = argv[3];
	
	/**
		* Read some meta data from the measurement set
		*/
	casa::MeasurementSet ms(msName);
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::MSAntenna aTable = ms.antenna();
	size_t antennaCount = aTable.nrow();
	if(antennaCount == 0) throw std::runtime_error("No antennae in set");
	casa::MPosition::ROScalarColumn antPosColumn(aTable, aTable.columnName(casa::MSAntennaEnums::POSITION));
	casa::MPosition ant1Pos = antPosColumn(0);
	casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	casa::MEpoch time = timeColumn(0);
	
	casa::MSField fTable(ms.field());
	if(fTable.nrow() != 1) throw std::runtime_error("Need exactly one field in set");
	casa::MDirection::ROScalarColumn refDirColumn(fTable, fTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	casa::MDirection refDir = refDirColumn(0);
	casa::MeasFrame frame(ant1Pos, time);
	const casa::MDirection::Ref hadecRef(casa::MDirection::HADEC, frame);
	const casa::MDirection::Ref azelgeoRef(casa::MDirection::AZELGEO, frame);
	const casa::MDirection::Ref j2000Ref(casa::MDirection::J2000, frame);
	casa::MPosition wgs = casa::MPosition::Convert(ant1Pos, casa::MPosition::WGS84)();
	double arrLatitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	casa::MDirection zenith(casa::MVDirection(0, M_PI), azelgeoRef);
	casa::MDirection zenithHaDec = casa::MDirection::Convert(zenith, hadecRef)();
	double zenithHa = zenithHaDec.getValue().getValue()[0];
	double zenithDec = zenithHaDec.getValue().getValue()[1];
	
	FitsReader reader(inpFitsname);
	size_t width = reader.ImageWidth();
	size_t height = reader.ImageHeight();
	double pixelSizeX = reader.PixelSizeX();
	double pixelSizeY = reader.PixelSizeY();
	
	std::vector<double> outImageX(width*height), outImageY(width*height);
	const casa::Unit radUnit("rad");
	
	casa::Table mwaTilePointing = ms.keywordSet().asTable("MWA_TILE_POINTING");
	casa::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
	casa::Array<int> delaysArr = delaysCol(0);
	casa::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
	double delays[16];
	std::cout << "Delays: [";
	for(int i=0; i!=16; ++i)
	{
		delays[i] = delaysArrPtr[i];
		std::cout << delays[i];
		if(i != 15) std::cout << ',';
	}
	std::cout << "]\n";
	//const double delays[16] = {6,9,12,15,4,7,10,13,2,5,8,11,0,3,6,9};
		
	TileBeam tilebeam(delays);
	
	double *xPtr = &outImageX[0], *yPtr = &outImageY[0];
	casa::MDirection::Convert
		j2000ToHaDecRef(j2000Ref, hadecRef),
		j2000ToAzelGeoRef(j2000Ref, azelgeoRef);
	for(size_t y=0;y!=height;++y)
	{
		for(size_t x=0;x!=width;++x)
		{
			double l, m, ra, dec;
			ImageCoordinates::XYToLM(x, y, pixelSizeX, pixelSizeY, width, height, l, m);
			ImageCoordinates::LMToRaDec(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
			
			double xPow, yPow;
			tilebeam.AnalyticGain(ra, dec, j2000Ref, j2000ToHaDecRef, j2000ToAzelGeoRef, arrLatitude, frequency, xPow, yPow);
			
			*xPtr = xPow;
			*yPtr = yPow;
			
			++xPtr; ++yPtr;
		}
		std::cout << '.' << std::flush;
	}
	
	std::cout << "\nWriting...\n";
	FitsWriter xWriter("beam-x.fits"), yWriter("beam-y.fits");
	
	xWriter.Write<double>(&outImageX[0], width, height, reader.PhaseCentreRA(), reader.PhaseCentreDec(), reader.PixelSizeX(), reader.PixelSizeY(), reader.Frequency(), reader.Bandwidth(), reader.DateObs());
	yWriter.Write<double>(&outImageY[0], width, height, reader.PhaseCentreRA(), reader.PhaseCentreDec(), reader.PixelSizeX(), reader.PixelSizeY(), reader.Frequency(), reader.Bandwidth(), reader.DateObs());
	
}
