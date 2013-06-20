#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ScalarColumn.h>

#include "sourcesdf.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: apparently [-ms <ms>] [-beam <beamfitsfile>] [-plot] <model> <fitsfile>\n"
			"Prints apparent model to stdout, constructed from given model and\n"
			"image fitsfile.\n"
			"-ms: read meta data from ms (e.g., time for -plot) will be used.\n"
			"-plot: don't write as model file, but as gnuplottable text file.\n"
			"-beam: read the beam and apply it to the found apparent values.\n";
	} else {
		int argi = 1;
		std::string msFilename, beamFilename;
		bool asPlot = false;
		while(argv[argi][0] == '-')
		{
			const char *option = &argv[argi][1];
			if(strcmp(option, "ms") == 0)
			{
				++argi;
				msFilename = argv[argi];
			}
			else if(strcmp(option, "beam") == 0)
			{
				++argi;
				beamFilename = argv[argi];
			}
			else if(strcmp(option, "plot") == 0)
			{
				asPlot = true;
			}
			else
			{
				throw std::runtime_error("Unable to parse options");
			}
			++argi;
		}
		
		const char *modelFilename = argv[argi];
		const char *fitsFilename = argv[argi+1];
		
		Model model(modelFilename);
		FitsReader fitsReader(fitsFilename);
		double *image = new double[fitsReader.ImageWidth() * fitsReader.ImageHeight()];
		fitsReader.Read(image);
		
		std::string setDescription;
		if(!msFilename.empty())
		{
			casa::MeasurementSet ms(msFilename);
			casa::ROScalarColumn<double> timeCol(ms, ms.columnName(casa::MeasurementSet::TIME));
			double startTime = timeCol(0);
			double endTime = timeCol(ms.nrow()-1);
			std::stringstream s;
			s << std::setprecision(12) << (endTime + startTime) * 0.5;
			setDescription = s.str();
		} else {
			setDescription = fitsFilename;
		}
		
		std::vector<double> beamData;
		if(!beamFilename.empty())
		{
			FitsReader beamReader(beamFilename);
			beamData.resize(beamReader.ImageWidth() * beamReader.ImageHeight());
			beamReader.Read(&beamData[0]);
		}

		if(asPlot)
		{
			std::cout << setDescription;
		}
		
		for(Model::const_iterator s=model.begin();s!=model.end();++s)
		{
			ModelSource source = *s;
			long double l, m;
			ImageCoordinates::RaDecToLM<long double>(source.PosRA(), source.PosDec(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), l, m);
			double
				x = l / fitsReader.PixelSizeX() + fitsReader.ImageWidth()/2,
				y = m / fitsReader.PixelSizeY() + fitsReader.ImageHeight()/2;
			double value = 0.0, beamValue = 1.0;
			//std::cout << x << ',' << y << '\n';
			if(x > 0.0 && y > 0.0 && x < fitsReader.ImageWidth() && y < fitsReader.ImageHeight())
			{
				size_t xi = size_t(floor(x)), yi = size_t(floor(y));
				value = image[yi*fitsReader.ImageWidth() + xi];
				
				++xi;
				if(xi < fitsReader.ImageWidth()) {
					value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
					
					++yi;
					if(yi < fitsReader.ImageHeight())
						value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
				}
				
				--xi;
				if(yi < fitsReader.ImageHeight())
					value = std::max(value, image[yi*fitsReader.ImageWidth() + xi]);
				
				if(!beamFilename.empty())
				{
					size_t xi = size_t(round(x)), yi = size_t(round(y));
					if(xi < fitsReader.ImageWidth() && yi < fitsReader.ImageHeight())
					{
						beamValue = beamData[yi*fitsReader.ImageWidth() + xi];
						//beamValue *= beamValue;
						value /= beamValue;
					}
				}
			}
			if(asPlot)
			{
				std::cout << '\t' << value;
				if(!beamFilename.empty())
				{
					std::cout << '\t' << beamValue;
				}
			}
			else {
				if(value > 0.0) {
					source.SetSED(SpectralEnergyDistribution(value, 1.0));
					std::cout << source.ToString() << '\n';
				}
			}
		}
		
		delete[] image;
		
		if(asPlot)
		{
			std::cout << '\n';
		}
	}
}
