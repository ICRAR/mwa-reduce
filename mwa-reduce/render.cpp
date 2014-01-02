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
		std::cout << "syntax: render [-delaunay] [-t templatefits] [-o <outputfits>] [-b] [-r] [-a] <model>\n";
	else {
		std::string templateFits;
		std::string outputFitsName;
		bool restore = false, addToTemplate = false, delaunay = false;
		
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
			else if(param == "delaunay") {
				delaunay = true;
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
		
		if(delaunay)
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
