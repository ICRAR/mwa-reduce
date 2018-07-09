#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include "model/measuredsed.h"
#include "model/model.h"

#include "fitsreader.h"

#include "units/imagecoordinates.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: apparently [-ms <ms>] [-beam <beamfitsfile>] [-plot] [-save <model>] <model> <fitsfile>\n"
			"Prints apparent model to stdout, constructed from given model and\n"
			"image fitsfile.\n"
			"-ms: read meta data from ms (e.g., time for -plot) will be used.\n"
			"-plot: don't write as model file, but as gnuplottable text file.\n"
			"-beam: read the beam and apply it to the found apparent values.\n"
			"-save: output as model.\n";
	} else {
		int argi = 1;
		std::string msFilename, beamFilename, outputFilename;
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
			else if(strcmp(option, "save") == 0)
			{
				++argi;
				outputFilename = argv[argi];
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
			casacore::MeasurementSet ms(msFilename);
			casacore::ROScalarColumn<double> timeCol(ms, ms.columnName(casacore::MeasurementSet::TIME));
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
		const size_t
			width = fitsReader.ImageWidth(),
			height = fitsReader.ImageHeight();
		
		Model measuredModel;
		for(Model::iterator s=model.begin();s!=model.end();++s)
		{
			ModelSource& source = *s;
			long double l, m;
			ImageCoordinates::RaDecToLM<long double>(source.Peak().PosRA(), source.Peak().PosDec(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), l, m);
			int x, y;
			ImageCoordinates::LMToXY<long double>(l, m, fitsReader.PixelSizeX(), fitsReader.PixelSizeY(), width, height, x, y);
			double value = 0.0, beamValue = 1.0;
			//std::cout << x << ',' << y << '\n';
			if(x >= 0 && y >= 0 && x < int(width) && y < int(height))
			{
				value = image[y*fitsReader.ImageWidth() + x];
				
				if(!beamFilename.empty())
				{
					if(x < int(width) && y < int(height))
					{
						beamValue = beamData[y*width + x];
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
					source.Peak().SetSED(MeasuredSED(value, fitsReader.Frequency()));
					measuredModel.AddSource(source);
				}
			}
		}
		
		delete[] image;
		
		if(asPlot)
		{
			std::cout << '\n';
		}
		
		if(!outputFilename.empty())
			measuredModel.Save(outputFilename);
	}
}
