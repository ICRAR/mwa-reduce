#include <iostream>
#include <stdexcept>
#include <vector>

#include "delaunay.h"
#include "fitsreader.h"
#include "model.h"
#include "modelrenderer.h"
#include "fitswriter.h"
#include "triangleinterpolator.h"
#include "imagecoordinates.h"

int main(int argc, char* argv[])
{
	if(argc == 1)
		std::cout << "syntax: render [-ion g/l/m <solutionfile>] [-t templatefits] [-o <outputfits>] [-b] [-r] [-a] <model>\n";
	else {
		std::string templateFits;
		std::string outputFitsName;
		std::string ionParameter, ionSolutionFilename;
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
				ionParameter = argv[argi];
				++argi;
				ionSolutionFilename = argv[argi];
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
		double ra = 0.0, dec = 0.0;
		double pixelSizeX = 0.01, pixelSizeY = 0.01;
		double bandwidth = 1000000.0, dateObs = 0.0, frequency = 150000000.0;
		double beamSize = 10.0*(M_PI/180.0/60.0);
		
		std::unique_ptr<FitsWriter> writer;
		std::vector<double> image;
		if(!templateFits.empty())
		{
			FitsReader reader(templateFits);
			width = reader.ImageWidth();
			height = reader.ImageHeight();
			image.resize(width * height);
			ra = reader.PhaseCentreRA();
			dec = reader.PhaseCentreDec();
			pixelSizeX = reader.PixelSizeX();
			pixelSizeY = reader.PixelSizeY();
			bandwidth = reader.Bandwidth();
			dateObs = reader.DateObs();
			frequency = reader.Frequency();
			if(reader.HasBeam())
				beamSize = reader.BeamMajorAxisRad();
			if(addToTemplate)
				reader.Read(&image[0]);
			
			writer.reset(new FitsWriter(reader));
		}
		
		if(!outputFitsName.empty())
		{
			ModelRenderer renderer(ra, dec, pixelSizeX, pixelSizeY);
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
			Delaunay triangulator;
			for(Model::iterator i=model.begin(); i!=model.end(); ++i)
			{
				ModelSource& source = *i;
				triangulator.AddVertex(source.MeanRA(), source.MeanDec(), &source);
			}
			triangulator.Triangulate();
			triangulator.SaveConvexHullAsKvis("convex.ann");
			triangulator.SaveTriangulationAsKvis("triangles.ann");
			
			TriangleInterpolator interpolator;
			for(size_t i=0; i!=triangulator.TriangleCount(); ++i)
			{
				Delaunay::Triangle triangle = triangulator.GetTriangle(i);
				ModelSource* source[3];
				double x[3], y[3];
				for(size_t j=0; j!=3; ++j) {
					source[j] = reinterpret_cast<ModelSource*>(triangle.userData[j]);
					double l, m;
					ImageCoordinates::RaDecToLM(triangle.x[j], triangle.y[j], ra, dec, l, m);
					ImageCoordinates::LMToXYfloat(l, m, pixelSizeX, pixelSizeY, width, height, x[j], y[j]);
					std::cout << x[j] << ',' << y[j] << " -> ";
				}
				std::cout << "Interpolating\n";
				interpolator.Interpolate(&image[0], width, height,
					x[0], y[0], source[0]->TotalFlux(frequency, 0),
					x[1], y[1], source[1]->TotalFlux(frequency, 0),
					x[2], y[2], source[2]->TotalFlux(frequency, 0)
				);
			}
			
			std::vector<double>
				xs(triangulator.ConvexVerticesCount()),
				ys(triangulator.ConvexVerticesCount());
			std::vector<ModelSource*>
				sources(triangulator.ConvexVerticesCount());
			for(size_t i=0; i!=triangulator.ConvexVerticesCount(); ++i)
			{
				double curra, curdec;
				triangulator.GetConvexVertex(i, curra, curdec, reinterpret_cast<void*&>(sources[i]));
				double l, m;
				ImageCoordinates::RaDecToLM(curra, curdec, ra, dec, l, m);
				ImageCoordinates::LMToXYfloat(l, m, pixelSizeX, pixelSizeY, width, height, xs[i], ys[i]);
			}
			for(size_t i=0; i!=triangulator.ConvexVerticesCount(); ++i)
			{
				//int i = 4;
				size_t
					a = i,
					b = (i+1) % triangulator.ConvexVerticesCount(),
					c = (i+2) % triangulator.ConvexVerticesCount();
				std::cout << "Interpolating " << xs[b] << ',' << ys[b] << '\n';
				interpolator.InterpolateConvexHullEdge(&image[0], width, height,
					xs[a], ys[a], sources[a]->TotalFlux(frequency, 0),
					xs[b], ys[b], sources[b]->TotalFlux(frequency, 0)
				);
				interpolator.InterpolateConvexHullVertex(&image[0], width, height,
					xs[a], ys[a],
					xs[b], ys[b], sources[b]->TotalFlux(frequency, 0),
					xs[c], ys[c]
				);
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
