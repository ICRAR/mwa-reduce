#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <set>

#include "sourcesdf.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "cleanalgorithms/simpleclean.h"

#include "areaset.h"
#include "parser/areaparser.h"
#include "uvector.h"

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: fitsmodel [options] <output model> <fitsfile> [spectral index] [ref freq MHz]\n"
			"Turns components in fitsfile into a model.\nOptions:\n"
			"\t-a <output areafile>\n"
			"\t-d <merge distance>\n"
			"\t-l <lower limit>\n"
			"\t-in <areafile>\n"
			"\t-out <areafile>\n";
	} else {
		size_t argi = 1;
		bool merge = false;
		double mergeDistance = 0.0;
		double limit = 0.0;
		std::string outputAreaFilename, insideFilename, outsideFilename;
		while(argv[argi][0] == '-')
		{
			std::string option(&argv[argi][1]);
			if(option == "a") {
				++argi;
				outputAreaFilename = argv[argi];
			}
			else if(option == "d") {
				++argi;
				merge = true;
				mergeDistance = atof(argv[argi]) * (M_PI / 180.0);
			}
			else if(option == "l") {
				++argi;
				limit = atof(argv[argi]);
			}
			else if(option == "in") {
				++argi;
				insideFilename = argv[argi];
			}
			else if(option == "out") {
				++argi;
				outsideFilename = argv[argi];
			}
			else throw std::runtime_error("Invalid param");
			++argi;
		}
		const char *modelFilename = argv[argi];
		++argi;
		const char *fitsFilename = argv[argi];
		long double spectralIndex, refFreq;
		if(argc - argi > 1)
			spectralIndex = atof(argv[argi + 1]);
		else
			spectralIndex = 0.0;
		FitsReader fitsReader(fitsFilename);
		if(argc - argi > 2)
			refFreq = atof(argv[argi + 2]) * 1000000;
		else
			refFreq = fitsReader.Frequency();
		const size_t
			width = fitsReader.ImageWidth(),
			height = fitsReader.ImageHeight();
		ao::uvector<double> image(width * height);
		fitsReader.Read(&image[0]);
		
		if(!insideFilename.empty() || !outsideFilename.empty())
		{
			std::cout << "Applying areas... " << std::flush;
			AreaParser parser;
			AreaSet areaIn, areaOut;
			if(!insideFilename.empty())
				parser.Parse(areaIn, insideFilename);
			if(!outsideFilename.empty())
				parser.Parse(areaOut, outsideFilename);
			for(size_t y=0; y!=height; ++y)
			{
				double* imageRow = &image[y*width];
				for(size_t x=0; x!=width; ++x)
				{
					if(imageRow[x] != 0.0)
					{
						long double l, m;
						ImageCoordinates::XYToLM<long double>(x, y, fitsReader.PixelSizeX(), fitsReader.PixelSizeY(), width, height, l, m);
						l += fitsReader.PhaseCentreDL(); m += fitsReader.PhaseCentreDM();
					
						ModelComponent component;
						long double ra, dec;
						ImageCoordinates::LMToRaDec<long double>(l, m, fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), ra, dec);
						
						bool isIn = insideFilename.empty() || areaIn.IsInArea(ra, dec);
						bool isOut = areaOut.IsInArea(ra, dec);
						if(!isIn || isOut)
							imageRow[x] = 0.0;
					}
				}
			}
			std::cout << "DONE\n";
		}
		
		Model model;
		CleanAlgorithm::GetModelFromImage(model, &image[0], width, height, fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), fitsReader.PixelSizeX(), fitsReader.PixelSizeY(), spectralIndex, refFreq);
		
		model.SortOnBrightness();
		
		if(merge)
		{
			std::set<size_t> sourcesToRemove;
			for(size_t i=0; i!=model.SourceCount(); ++i)
			{
				ModelSource &refSource = model.Source(i);
				for(size_t j=i+1; j!=model.SourceCount(); ++j)
				{
					const ModelSource &source2 = model.Source(j);
					double distance = 
						ImageCoordinates::AngularDistance(refSource.Peak().PosRA(), refSource.Peak().PosDec(), source2.Peak().PosRA(), source2.Peak().PosDec());
					if(distance <= mergeDistance)
					{
						sourcesToRemove.insert(j);
					}
				}
			}
			for(std::set<size_t>::reverse_iterator i=sourcesToRemove.rbegin();
					i!=sourcesToRemove.rend(); ++i)
			{
				model.RemoveSource(*i);
			}
		}
		
		if(limit > 0.0)
		{
			std::cout << "Thresholding sources... " << std::flush;
			double totalFlux = 0.0, removedFlux = 0.0;
			for(size_t i=model.SourceCount(); i>0; --i)
			{
				size_t sIndex = i - 1;
				double flux = model.Source(sIndex).TotalFlux(refFreq, Polarization::StokesI);
				if(flux < limit)
				{
					model.RemoveSource(sIndex);
					removedFlux += flux;
				}
				totalFlux += flux;
			}
			std::cout << "DONE (removed: " << removedFlux << " / " << totalFlux << " Jy)\n";
		}
		
		AreaSet areaSet;
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			ModelSource &source = model.Source(i);
			std::ostringstream str;
			str << "component" << (i+1);
			source.SetName(str.str());
			
			SkyArea area;
			SkyAreaElement element;
			element.SetToCircle(mergeDistance, source.Peak().PosRA(), source.Peak().PosDec());
			area.AddElement(element);
			area.SetName(source.Name());
			areaSet.AddArea(area);
		}
		if(!outputAreaFilename.empty())
		{
			std::ofstream areaFile(outputAreaFilename.c_str());
			areaSet.Save(areaFile);
		}
		std::cout << "Writing model with " << model.SourceCount() << " sources... " << std::flush;
		model.Save(modelFilename);
		std::cout << "DONE\n";
	}
}
