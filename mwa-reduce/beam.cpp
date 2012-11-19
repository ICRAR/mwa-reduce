#include "tilebeam.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "fitswriter.h"

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

int main(int argc, char *argv[])
{
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
	double latitude = wgs.getValue().getLat(); // ant1Pos.getValue().getLat();
	
	FitsReader reader(inpFitsname);
	size_t width = reader.ImageWidth();
	size_t height = reader.ImageHeight();
	double pixelSizeX = reader.PixelSizeX();
	double pixelSizeY = reader.PixelSizeY();
	
	std::vector<double> outImageX(width*height), outImageY(width*height);
	const casa::Unit radUnit("rad");
	double midX = (double) width / 2.0, midY = (double) height / 2.0;
	const double delays[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	TileBeam tilebeam(delays);
	
	double *xPtr = &outImageX[0], *yPtr = &outImageY[0];
	casa::MDirection imageDir(casa::MVDirection(
		casa::Quantity(0.0, radUnit),     // RA
		casa::Quantity(0.0, radUnit)),  // DEC
		j2000Ref);
	casa::MDirection::Convert
		j2000ToHaDec(j2000Ref, hadecRef),
		j2000ToAzelGeo(j2000Ref, azelgeoRef);
	for(size_t y=0;y!=height;++y)
	{
		for(size_t x=0;x!=width;++x)
		{
			double l = ((double) x - midX) * pixelSizeX;
			double m = (midY - (double) y) * pixelSizeY;
			double ra, dec;
			ImageCoordinates::LMToRaDec(l, m, reader.PhaseCentreRA(), reader.PhaseCentreDec(), ra, dec);
			
			imageDir.set(casa::MVDirection(
				casa::Quantity(ra, radUnit),     // RA
				casa::Quantity(dec,radUnit)));
			// convert ra, dec to za, az
			casa::MDirection hadec = j2000ToHaDec(imageDir);
			double ha = hadec.getValue().get()[0];
			double sinLat, cosLat;
			sincos(latitude, &sinLat, &cosLat);
			double sinDec, cosDec;
			sincos(dec, &sinDec, &cosDec);
			double cosHA = cos(ha);
			double zenithDistance = acos(sinLat * sinDec + cosLat * cosDec * cosHA);
			
			casa::MDirection azel = j2000ToAzelGeo(imageDir);
			double azimuth = azel.getValue().get()[0];
			
			double xPow, yPow;
			tilebeam.AnalyticGain(zenithDistance, azimuth, frequency, xPow, yPow);
			
			*xPtr = xPow;
			*yPtr = yPow;
			
			++xPtr; ++yPtr;
		}
		std::cout << '.' << std::flush;
	}
	
	FitsWriter xWriter("beam-x.fits"), yWriter("beam-y.fits");
	
	xWriter.Write<double>(&outImageX[0], width, reader.PhaseCentreRA(), reader.PhaseCentreDec(), reader.PixelSizeX(), reader.PixelSizeY());
	yWriter.Write<double>(&outImageY[0], width, reader.PhaseCentreRA(), reader.PhaseCentreDec(), reader.PixelSizeX(), reader.PixelSizeY());
}
