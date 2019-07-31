#include <iostream>
#include <stdexcept>
#include <vector>
#include <random>

#include "fitsreader.h"
#include "fitswriter.h"
#include "ioninterpolator.h"
#include "image.h"
#include "modelrenderer.h"
#include "banddata.h"

#include "units/angle.h"

#include "model/model.h"

void meanPos(const std::vector<ModelSource*>& sources, double& ra, double& dec);

int main(int argc, char* argv[])
{
	if(argc == 1)
		std::cout << "syntax: render [-n <noiselevel>] [-ion <solutionfile> <outprefix>] [-t templatefits] [-o <outputfits>] [-b] [-r [-beam <maj> <min> <pa>]] [-a] [-centre <ra> <dec>] [-size <width> <height>] [-scale <scale>] [-frequency <valueHz>] <model>\n";
	else {
		std::string templateFits;
		std::string outputFitsName;
		std::string ionOutPrefix;
		const char* ionSolutionFilename = 0;
		bool restore = false, addToTemplate = false, ionospheric = false, hasManualBeam = false;
		int argi = 1;
		double ra = 0.0, dec = 0.0, dl = 0.0, dm = 0.0;
		double pixelSizeX = 0.012*(M_PI/180.0), pixelSizeY = 0.012*(M_PI/180.0);
		double noise = 0.0;
		double
			beamMaj = 2.0*(M_PI/180.0/60.0),
			beamMin = 2.0*(M_PI/180.0/60.0),
			beamPA = 0.0;
		size_t sizeWidth = 0, sizeHeight = 0;
		double setFrequency = 0.0;
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
			else if(param == "n") {
				++argi;
				noise = atof(argv[argi]);
			}
			else if(param == "beam") {
				hasManualBeam = true;
				++argi;
				beamMaj = Angle::Parse(argv[argi], "beam major axis", Angle::Arcseconds);
				++argi;
				beamMin = Angle::Parse(argv[argi], "beam minor axis", Angle::Arcseconds);
				++argi;
				beamPA = Angle::Parse(argv[argi], "beam position angle", Angle::Degrees);
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
			else if(param == "centre")
			{
				++argi;
				ra = RaDecCoord::ParseRA(argv[argi]);
				++argi;
				dec = RaDecCoord::ParseDec(argv[argi]);
			}
			else if(param == "size")
			{
				sizeWidth = atoi(argv[argi+1]);
				sizeHeight = atoi(argv[argi+2]);
				argi += 2;
			}
			else if(param == "scale")
			{
				++argi;
				pixelSizeX = Angle::Parse(argv[argi], "scale", Angle::Degrees);
				pixelSizeY = pixelSizeX;
			}
			else if(param == "frequency")
			{
				++argi;
				setFrequency = atof(argv[argi]);
			}
			else throw std::runtime_error("Invalid param");
			++argi;
		}
	
		Model model(argv[argi]);
	
		size_t width = 4096, height = 4096;
		double bandwidth = 1000000.0, dateObs = 0.0, frequency = 150000000.0;
		
		ImageBufferAllocator allocator;
		std::unique_ptr<FitsWriter> writer;
		std::unique_ptr<FitsReader> reader;
		Image image;
		if(!templateFits.empty())
		{
			double wscImgWeight = 0.0;
			reader.reset(new FitsReader(templateFits));
			width = reader->ImageWidth();
			height = reader->ImageHeight();
			image = Image(width, height, allocator);
			ra = reader->PhaseCentreRA();
			dec = reader->PhaseCentreDec();
			dl = reader->PhaseCentreDL();
			dm = reader->PhaseCentreDM();
			pixelSizeX = reader->PixelSizeX();
			pixelSizeY = reader->PixelSizeY();
			bandwidth = reader->Bandwidth();
			dateObs = reader->DateObs();
			frequency = reader->Frequency();
			if(reader->HasBeam() && !hasManualBeam)
			{
				beamMaj = reader->BeamMajorAxisRad();
				beamMin = reader->BeamMinorAxisRad();
				beamPA = reader->BeamPositionAngle();
			}
			reader->ReadDoubleKeyIfExists("WSCIMGWG", wscImgWeight);
			if(addToTemplate)
				reader->Read(&image[0]);
			
			writer.reset(new FitsWriter(*reader));
			if(wscImgWeight != 0.0)
				writer->SetExtraKeyword("WSCIMGWG", wscImgWeight);
		}
		else {
			image = Image(width, height, allocator);
			writer.reset(new FitsWriter());
		}
		
		if(sizeWidth!=0 && sizeHeight!=0)
		{
			if(sizeWidth > width && sizeHeight > height)
			{
				image = image.Untrim(sizeWidth, sizeHeight);
				width = sizeWidth;
				height = sizeHeight;
			}
		}
		
		if(setFrequency != 0.0)
		{
			frequency = setFrequency;
			writer->SetFrequency(setFrequency, writer->Bandwidth());
		}
			
		if(!outputFitsName.empty())
		{
			ModelRenderer renderer(ra, dec, pixelSizeX, pixelSizeY, dl, dm);
			if(noise != 0.0)
			{
				std::random_device rd;
				std::mt19937 rnd(rd());
				std::normal_distribution<double> dist(0.0, noise);
				for(size_t i=0; i!=width*height; ++i)
					image[i] += dist(rnd);
			}
			if(restore)
			{
				renderer.Restore(&image[0], width, height, model, beamMaj, beamMin, beamPA, frequency-bandwidth*0.5, frequency+bandwidth*0.5, Polarization::StokesI);
			}
			else {
				renderer.RenderModel(&image[0], width, height, model, frequency-bandwidth*0.5, frequency+bandwidth*0.5, Polarization::StokesI);
			}
		}
		
		if(ionospheric)
		{
			IonInterpolator interpolator(model, ra, dec, pixelSizeX, pixelSizeY, width, height);
			IonSolutionFile solutionFile;
			solutionFile.OpenForReading(ionSolutionFilename);
			std::cout << "Model has " << model.ComponentCount() << " components in " << model.ClusterCount() << " clusters.\n";
			std::cout << "Directions in solution file: " << solutionFile.DirectionCount() << '\n';
			if(solutionFile.DirectionCount() != model.ClusterCount())
				throw std::runtime_error("Direction count in solutions and cluster count in model do not match!");
			interpolator.InitializeFile(solutionFile);
			
			std::ofstream
				plotPosFile(ionOutPrefix+"-posplot.plt"),
				plotGainFile(ionOutPrefix+"-gainplot.plt"),
				vectorPlot(ionOutPrefix+"-vectors.ann");
			plotPosFile
				<< "set terminal postscript enhanced color\n"
				<< "#set logscale y\n"
				<< "#set xrange [0.001:]\n"
				<< "#set yrange [-8:2]\n"
				<< "set output \"" << ionOutPrefix << "-posplot.ps\"\n"
				<< "set key bottom left\n"
				<< "set xlabel \"Freq (\\lambda^2)\"\n"
				<< "set ylabel \"Pos offset (arcmin)\"\n"
				<< "plot\\\n";
			plotGainFile
				<< "set terminal postscript enhanced color\n"
				<< "#set logscale y\n"
				<< "#set xrange [0.001:]\n"
				<< "#set yrange [-8:2]\n"
				<< "set output \"" << ionOutPrefix << "-gainplot.ps\"\n"
				<< "set key bottom left\n"
				<< "set xlabel \"Freq (\\lambda^2)\"\n"
				<< "set ylabel \"Gain\"\n"
				<< "plot\\\n";
				
			for(size_t s=0; s!=solutionFile.DirectionCount(); ++s)
			{
				double meanRA, meanDec;
				interpolator.GetMeanPosForDirection(s, meanRA, meanDec);
				
				std::ostringstream name;
				name << ionOutPrefix << "-direc" << s << "-";
				if(s != 0) {
					plotPosFile << ",\\\n";
					plotGainFile << ",\\\n";
				}
				std::ofstream
					direcfile(name.str() + "posoff.txt"),
					gainfile(name.str() + "gain.txt");
				plotPosFile << "\"" << name.str() << "posoff.txt\" using 2:3 with points title \"\"";
				plotGainFile << "\"" << name.str() << "gain.txt\" using 2:3 with points title \"\"";
				double totalDl = 0.0, totalDm = 0.0;
				size_t totalCount = 0;
				for(size_t cb=0; cb!=solutionFile.ChannelBlockCount(); ++cb)
				{
					double freq = (solutionFile.EndFrequency() - solutionFile.StartFrequency()) *
						double(cb) / solutionFile.ChannelBlockCount() + solutionFile.StartFrequency();
					double lambda = BandData::FrequencyToLambda(freq);
					double sumDl = 0.0, sumDm = 0.0, sumG = 0.0;
					size_t count = 0;
					for(size_t interval=0; interval!=solutionFile.IntervalCount(); ++interval)
					{
						double
							dl = solutionFile.ReadSolution(IonSolutionFile::DlSolution, interval, cb, 0, s),
							dm = solutionFile.ReadSolution(IonSolutionFile::DmSolution, interval, cb, 0, s),
							g = solutionFile.ReadSolution(IonSolutionFile::GainSolution, interval, cb, 0, s);
						if(std::isfinite(dl) && std::isfinite(dm) && std::isfinite(g))
						{
							sumDl += dl; sumDm += dm; sumG += g;
							++count;
						}
					}
					totalDl += sumDl;
					totalDm += sumDm;
					totalCount += count;
					double
						posOff = sqrt(sumDl*sumDl + sumDm*sumDm) / count,
						gain = sumG / count;
					posOff *= 60.0*180.0/M_PI; // to arcmin
					direcfile << freq << '\t' << lambda*lambda << '\t' << posOff << '\n';
					gainfile << freq << '\t' << lambda*lambda << '\t' << gain << '\n';
				}
				totalDl *= 100.0/totalCount;
				totalDm *= 100.0/totalCount;
				vectorPlot
					<< "LINE\t"
					<< meanRA*(180.0/M_PI) << '\t' << meanDec*(180.0/M_PI) << '\t'
					<< (meanRA-totalDl)*(180.0/M_PI) << '\t'
					<< (meanDec-totalDm)*(180.0/M_PI) << '\n';
			}
			plotPosFile << '\n';
			plotGainFile << '\n';
			vectorPlot.close();
			
			ao::uvector<double> interpolatedImage(width * height);
			for(size_t interval=0; interval!=solutionFile.IntervalCount(); ++interval)
			{
				std::cout << "Rendering interval " << interval << "...\n";
				interpolator.Initialize(solutionFile, interval, interval+1, 0, solutionFile.ChannelBlockCount(), 0);
				
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
