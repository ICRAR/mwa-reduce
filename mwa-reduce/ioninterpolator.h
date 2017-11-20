#ifndef ION_INTERPOLATOR_H
#define ION_INTERPOLATOR_H

#include "delaunay.h"
#include "units/imagecoordinates.h"
#include "ionsolutionfile.h"
#include "model/model.h"
#include "triangleinterpolator.h"
#include "fitsreader.h"
#include "uvector.h"

class IonInterpolator
{
private:
	typedef std::size_t size_t;
public:
	IonInterpolator(const Model& model, const FitsReader& templateReader) :
		_model(model),
		_directionIndices(),
		_phaseCentreRA(templateReader.PhaseCentreRA()),
		_phaseCentreDec(templateReader.PhaseCentreDec()),
		_pixelSizeX(templateReader.PixelSizeX()),
		_pixelSizeY(templateReader.PixelSizeY()),
		_width(templateReader.ImageWidth()),
		_height(templateReader.ImageHeight())
	{
	}
	
	IonInterpolator(const Model& model, double ra, double dec, double scaleX, double scaleY, size_t width, size_t height) :
		_model(model),
		_directionIndices(),
		_phaseCentreRA(ra),
		_phaseCentreDec(dec),
		_pixelSizeX(scaleX),
		_pixelSizeY(scaleY),
		_width(width),
		_height(height)
	{
	}
	
	void InitializeFile(IonSolutionFile& solutions)
	{
		_directionIndices.resize(solutions.DirectionCount());
		for(size_t i=0; i!=solutions.DirectionCount(); ++i)
			_directionIndices[i] = i;
		solutions.ReadClusterMetaInfo(_model, _sourcesPerDirection);
	}
	
	void GetMeanPosForDirection(size_t direction, double& meanRA, double& meanDec)
	{
		meanPos(_sourcesPerDirection[direction], meanRA, meanDec);
	}
	
	/**
	 * Initialize the ioninterpolator to be able to interpolate solutions
	 * for a given ionospheric interval, channel and polarization.
	 */
	void Initialize(IonSolutionFile& solutions, size_t startInterval, size_t endInterval, size_t startChannel, size_t endChannel, size_t polarization)
	{
		InitializeFile(solutions);
	
		_startInterval = startInterval;
		_endInterval = endInterval;
		_startChannel = startChannel;
		_endChannel = endChannel;
		_polarization = polarization;
		
		initializeTriangulator(_triangulator, solutions);
		
	}
	
	void initializeTriangulator(Delaunay& delaunay, IonSolutionFile& solutions)
	{
		delaunay.Clear();
		for(size_t direction=0; direction!=solutions.DirectionCount(); ++direction)
		{
			double posRA, posDec;
			GetMeanPosForDirection(direction, posRA, posDec);
			if(!std::isfinite(posRA) || !std::isfinite(posDec))
				throw std::runtime_error("One of the directions of the solutions has a non-finite RA or Dec");
			IonSolutionFile::Solution solution;
			getAverageSolution(solutions, solution, direction);
			if(!std::isfinite(solution.dl) || !std::isfinite(solution.dm) || !std::isfinite(solution.gain))
			{
				const ModelSource* source = _sourcesPerDirection[direction][0];
				std::cout << "Warning: direction " << source->ClusterName() << " did not have finite average solution over interval " << _startInterval << " to " << _endInterval << ", channel " << _startChannel << " to " << _endChannel << "!\n";
			}
			// Invert the ionospheric corrections: the source in our model
			//   is observed at position l+dl, m+dm, so if we want to know
			//   what is really at position l,m, we need to displace
			//   the observed source at l+dl, m+dm to l,m. Therefore, at position l+dl,m+dm,
			//   the ionospheric correction is dl,dm.
			double l, m;
			ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, l, m);
			l += solution.dl;
			m += solution.dm;
			ImageCoordinates::LMToRaDec(l, m, _phaseCentreRA, _phaseCentreDec, posRA, posDec);
			if(!std::isfinite(posRA) || !std::isfinite(posDec))
				throw std::runtime_error("One of the directions in the model is not finite.");
			delaunay.AddVertex(posRA, posDec, &_directionIndices[direction]);
		}
		delaunay.Triangulate();
		delaunay.SaveConvexHullAsKvis("convex.ann");
		delaunay.SaveTriangulationAsKvis("triangles.ann");
	}
	
	
	void Interpolate(double* image, IonSolutionFile& file, IonSolutionFile::IonSolutionType solutionType)
	{
		ao::uvector<double> solutionValues(file.DirectionCount());
		for(size_t direction=0; direction!=file.DirectionCount(); ++direction)
		{
			solutionValues[direction] = getAverageSolution(file, solutionType, direction);
		}
		
		TriangleInterpolator interpolator;
		for(size_t i=0; i!=_triangulator.TriangleCount(); ++i)
		{
			Delaunay::Triangle triangle = _triangulator.GetTriangle(i);
			size_t indices[3];
			double x[3], y[3];
			for(size_t j=0; j!=3; ++j) {
				indices[j] = *reinterpret_cast<size_t*>(triangle.userData[j]);
				double l, m;
				ImageCoordinates::RaDecToLM(triangle.x[j], triangle.y[j], _phaseCentreRA, _phaseCentreDec, l, m);
				ImageCoordinates::LMToXYfloat(l, m, _pixelSizeX, _pixelSizeY, _width, _height, x[j], y[j]);
				std::cout << x[j] << ',' << y[j] << " -> ";
			}
			std::cout << "Interpolating\n";
			interpolator.Interpolate(&image[0], _width, _height,
				x[0], y[0], solutionValues[indices[0]],
				x[1], y[1], solutionValues[indices[1]],
				x[2], y[2], solutionValues[indices[2]]
			);
		}
		
		std::vector<double>
			xs(_triangulator.ConvexVerticesCount()),
			ys(_triangulator.ConvexVerticesCount());
		std::vector<size_t> indices(_triangulator.ConvexVerticesCount());
		for(size_t i=0; i!=_triangulator.ConvexVerticesCount(); ++i)
		{
			Delaunay::ConvexVertex v = _triangulator.GetConvexVertex(i);
			indices[i] = *reinterpret_cast<size_t*>(v.userData);
			double l, m;
			ImageCoordinates::RaDecToLM(v.x, v.y, _phaseCentreRA, _phaseCentreDec, l, m);
			ImageCoordinates::LMToXYfloat(l, m, _pixelSizeX, _pixelSizeY, _width, _height, xs[i], ys[i]);
		}
		for(size_t i=0; i!=_triangulator.ConvexVerticesCount(); ++i)
		{
			//int i = 4;
			size_t
				a = i,
				b = (i+1) % _triangulator.ConvexVerticesCount(),
				c = (i+2) % _triangulator.ConvexVerticesCount();
			std::cout << "Interpolating " << xs[b] << ',' << ys[b] << '\n';
			interpolator.InterpolateConvexHullEdge(&image[0], _width, _height,
				xs[a], ys[a], solutionValues[indices[a]],
				xs[b], ys[b], solutionValues[indices[b]]
			);
			interpolator.InterpolateConvexHullVertex(&image[0], _width, _height,
				xs[a], ys[a],
				xs[b], ys[b], solutionValues[indices[b]],
				xs[c], ys[c]
			);
		}
	}
private:
	void getAverageSolution(IonSolutionFile& solutionFile, IonSolutionFile::Solution& solution, size_t direction) const
	{
		solution.dl = 0.0;
		solution.dm = 0.0;
		solution.gain = 0.0;
		size_t count = 0;
		for(size_t interval=_startInterval; interval!=_endInterval; ++interval)
		{
			for(size_t ch=_startChannel; ch!=_endChannel; ++ch)
			{
				IonSolutionFile::Solution s;
				solutionFile.ReadSolution(s, interval, ch, _polarization, direction);
				if(std::isfinite(s.dl) && std::isfinite(s.dm) && std::isfinite(s.gain))
				{
					solution.dl += s.dl;
					solution.dm += s.dm;
					solution.gain += s.gain;
					++count;
				}
			}
		}
		solution.dl /= count;
		solution.dm /= count;
		solution.gain /= count;
	}
	
	double getAverageSolution(IonSolutionFile& solutionFile, IonSolutionFile::IonSolutionType stype, size_t direction) const
	{
		size_t count = 0;
		double val = 0.0;
		for(size_t interval=_startInterval; interval!=_endInterval; ++interval)
		{
			for(size_t ch=_startChannel; ch!=_endChannel; ++ch)
			{
				double s = solutionFile.ReadSolution(stype, interval, ch, _polarization, direction);
				if(std::isfinite(s))
				{
					val += s;
					++count;
				}
			}
		}
		return val / count;
	}
	
	static void meanPos(const std::vector<ModelSource*>& sources, double& meanRA, double& meanDec)
	{
		std::vector<double> raValues;
		double decSum = 0.0;
		
		raValues.reserve(sources.size());
		for(std::vector<ModelSource*>::const_iterator s=sources.begin(); s!=sources.end(); ++s)
		{
			raValues.push_back((*s)->MeanRA());
			decSum += (*s)->MeanDec();
		}
		meanDec = decSum / sources.size();
		meanRA = ImageCoordinates::MeanRA(raValues);
	}

	Model _model;
	Delaunay _triangulator;
	std::vector<size_t> _directionIndices;

	const double
		_phaseCentreRA, _phaseCentreDec,
		_pixelSizeX, _pixelSizeY;
	const size_t
		_width, _height;
	
	size_t _startInterval, _endInterval, _startChannel, _endChannel, _polarization;
	std::vector<std::vector<ModelSource*>> _sourcesPerDirection;
};

#endif
