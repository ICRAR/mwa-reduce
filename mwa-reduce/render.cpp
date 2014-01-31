#include <iostream>
#include <stdexcept>
#include <vector>

#include "fitsreader.h"
#include "fitswriter.h"
#include "ioninterpolator.h"
#include "model.h"
#include "modelrenderer.h"

int main(int argc, char* argv[])
{
	if(argc == 1)
		std::cout << "syntax: render [-ion <solutionfile> <outprefix>] [-t templatefits] [-o <outputfits>] [-b] [-r] [-a] <model>\n";
	else {
		std::string templateFits;
		std::string outputFitsName;
		std::string ionOutPrefix;
		const char* ionSolutionFilename = 0;
		bool restore = false, addToTemplate = false, ionospheric = false;
		int argi = 1;
		while(argi < argc && argv[argi][0] == '-')
		{
			std::string param(&argv[argi][1]);
			if(param == "t") {
				++argi;
				templateFits = argv[argi];
			}
			else if(param == "r") {
				restore = true;
			}
			else if(param == "a") {
				addToTemplate = true;
			}
			else if(param == "ion") {
				ionospheric = true;
				++argi;
				ionSolutionFilename = argv[argi];
				++argi;
				ionOutPrefix = argv[argi];
			}
			else if(param == "o")
			{
				++argi;
				outputFitsName = argv[argi];
			}
			else throw std::runtime_error("Invalid param");
			++argi;
		}
	
		Model model(argv[argi]);
	
		size_t width = 1024, height = 1024;
		double ra = 0.0, dec = 0.0, dl = 0.0, dm = 0.0;
		double pixelSizeX = 0.01, pixelSizeY = 0.01;
		double bandwidth = 1000000.0, dateObs = 0.0, frequency = 150000000.0;
		double beamSize = 10.0*(M_PI/180.0/60.0);
		
		std::unique_ptr<FitsWriter> writer;
		std::unique_ptr<FitsReader> reader;
		std::vector<double> image;
		if(!templateFits.empty())
		{
			reader.reset(new FitsReader(templateFits));
			width = reader->ImageWidth();
			height = reader->ImageHeight();
			image.resize(width * height);
			ra = reader->PhaseCentreRA();
			dec = reader->PhaseCentreDec();
			dl = reader->PhaseCentreDL();
			dm = reader->PhaseCentreDM();
			pixelSizeX = reader->PixelSizeX();
			pixelSizeY = reader->PixelSizeY();
			bandwidth = reader->Bandwidth();
			dateObs = reader->DateObs();
			frequency = reader->Frequency();
			if(reader->HasBeam())
				beamSize = reader->BeamMajorAxisRad();
			if(addToTemplate)
				reader->Read(&image[0]);
			
			writer.reset(new FitsWriter(*reader));
		}
		
		if(!outputFitsName.empty())
		{
			ModelRenderer renderer(ra, dec, pixelSizeX, pixelSizeY, dl, dm);
			if(restore)
			{
				renderer.Restore(&image[0], width, height, model, beamSize, frequency-bandwidth*0.5, frequency+bandwidth*0.5, 0);
			}
			else {
				renderer.RenderModel(&image[0], width, height, model, frequency-bandwidth*0.5, frequency+bandwidth*0.5, 0);
			}
		}
		
		if(ionospheric)
		{
			IonInterpolator interpolator(model, *reader);
			IonSolutionFile solutionFile;
			solutionFile.OpenForReading(ionSolutionFilename);
			ao::uvector<double> interpolatedImage(width * height);
			for(size_t interval=0; interval!=solutionFile.IntervalCount(); ++interval)
			{
				std::cout << "Rendering interval " << interval << "...\n";
				interpolator.Initialize(solutionFile, interval, interval+1, 0, solutionFile.ChannelCount(), 0);
				
				std::ostringstream extStr;
				extStr << "-i";
				if(interval < 100)
				{
					extStr << '0';
					if(interval < 10) extStr << '0';
				}
				extStr << interval;
				extStr << ".fits";
				
				std::string gainfile = ionOutPrefix + "-gain" + extStr.str();
				interpolator.Interpolate(interpolatedImage.data(), solutionFile, IonSolutionFile::GainSolution);
				writer->Write(gainfile, interpolatedImage.data());
				
				std::string dlfile = ionOutPrefix + "-dl" + extStr.str();
				interpolator.Interpolate(interpolatedImage.data(), solutionFile, IonSolutionFile::DlSolution);
				writer->Write(dlfile, interpolatedImage.data());
				
				std::string dmfile = ionOutPrefix + "-dm" + extStr.str();
				interpolator.Interpolate(interpolatedImage.data(), solutionFile, IonSolutionFile::DmSolution);
				writer->Write(dmfile, interpolatedImage.data());
			}
		}
		
		if(!outputFitsName.empty())
		{
			writer->SetImageDimensions(width, height, ra, dec, pixelSizeX, pixelSizeY);
			writer->SetFrequency(frequency, bandwidth);
			writer->SetDate(dateObs);
			writer->Write(outputFitsName, &image[0]);
		}
	}
}
