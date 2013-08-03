#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <set>

#include "sourcesdf.h"
#include "model.h"
#include "fitsreader.h"
#include "imagecoordinates.h"
#include "cleanalgorithm.h"
#include "areaset.h"

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cout << "Usage: fitsmodel [-m <merge distance>] [-l <lower limit>] <fitsfile> [spectral index] [ref freq]\n"
			"Prints components in fitsfile to stdout as model.\n";
	} else {
		size_t argi = 1;
		bool merge = false;
		double mergeDistance = 0.0;
		double limit = 0.0;
		while(argv[argi][0] == '-')
		{
			std::string option(&argv[argi][1]);
			if(option == "m")
			{
				++argi;
				merge = true;
				mergeDistance = atof(argv[argi]) * (M_PI / 180.0);
			}
			else if(option == "l")
			{
				++argi;
				limit = atof(argv[argi]);
			}
			else throw std::runtime_error("Invalid param");
			++argi;
		}
		const char *fitsFilename = argv[argi];
		long double spectralIndex, refFreq;
		if(argc - argi > 1)
			spectralIndex = atof(argv[argi + 1]);
		else
			spectralIndex = 0.0;
		if(argc - argi > 2)
			refFreq = atof(argv[argi + 2]) * 1000000;
		else
			refFreq = 100000000.0;
		FitsReader fitsReader(fitsFilename);
		double *image = new double[fitsReader.ImageWidth() * fitsReader.ImageHeight()];
		fitsReader.Read(image);
		
		Model model;
		CleanAlgorithm::GetModelFromImage(model, image, fitsReader.ImageWidth(), fitsReader.ImageHeight(), fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(), fitsReader.PixelSizeX(), fitsReader.PixelSizeY(), spectralIndex, refFreq);
		
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
						ImageCoordinates::AngularDistance(refSource.PosRA(), refSource.PosDec(), source2.PosRA(), source2.PosDec());
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
			for(size_t i=model.SourceCount(); i>0; --i)
			{
				size_t sIndex = model.SourceCount() - 1;
				if(model.Source(sIndex).SED().FluxAtLowestFrequency() < limit)
					model.RemoveSource(i);
			}
		}
		
		AreaSet areaSet;
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			ModelSource &source = model.Source(i);
			std::ostringstream str;
			str << "component" << (i+1);
			source.SetName(str.str());
			std::cout << source.ToString() << '\n';
			
			SkyArea area;
			SkyAreaElement element;
			element.SetToCircle(mergeDistance, source.PosRA(), source.PosDec());
			area.AddElement(element);
			area.SetName(source.Name());
			areaSet.AddArea(area);
		}
		std::ofstream areaFile("model-areas.txt");
		areaSet.Save(areaFile);
		
		delete[] image;
	}
}
