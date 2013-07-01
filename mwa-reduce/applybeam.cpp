#include "fitsreader.h"
#include "fitswriter.h"

#include <vector>
#include <stdexcept>
#include <string>

int main(int argc, char *argv[])
{
	const char *inpFits = argv[1];
	const char *beamFits = argv[2];
	const char *outFits = argv[3];
	
	FitsReader inpReader(inpFits);
	size_t
		width = inpReader.ImageWidth(),
		height = inpReader.ImageHeight();
	FitsReader beamReader(beamFits);
	if(beamReader.ImageWidth() != width || beamReader.ImageHeight() != height)
		throw std::runtime_error("Beam and image do not have same size!");
	
	std::vector<double> inpImage(width*height), beamImage(width*height);
	
	inpReader.Read<double>(&inpImage[0]);
	beamReader.Read<double>(&beamImage[0]);
	
	std::vector<double>::iterator beamIter = beamImage.begin();
	for(std::vector<double>::iterator i=inpImage.begin(); i!=inpImage.end(); ++i)
	{
		*i /= *beamIter * *beamIter;
		++beamIter;
	}
	
	FitsWriter writer(outFits);
	writer.Write<double>(&inpImage[0], width, height, inpReader.PhaseCentreRA(), inpReader.PhaseCentreDec(), inpReader.PixelSizeX(), inpReader.PixelSizeY(), inpReader.Frequency(), inpReader.Bandwidth(), inpReader.DateObs());
}
