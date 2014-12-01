#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include "matrix2x2.h"

#include <cmath>
#include <cstring>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

#include "uvector.h"

class GaussianFitter
{
public:
	void Fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA)
	{
		size_t prefSize = std::max<size_t>(10, std::ceil(beamEst*10.0));
		if(prefSize%2 != 0) ++prefSize;
		if(prefSize < width || prefSize < height)
		{
			size_t boxWidth  = std::min(prefSize, width);
			size_t boxHeight = std::min(prefSize, height);
			size_t nIter = 0;
			bool boxWasLargeEnough;
			do {
				fit2DGaussianCentredInBox(image, width, height, beamEst, beamMaj, beamMin, beamPA, boxWidth, boxHeight);
				
				boxWasLargeEnough =
					(beamMaj*4.0 < boxWidth || width>=boxWidth) &&
					(beamMaj*4.0 < boxHeight || height>=boxHeight);
				if(!boxWasLargeEnough)
				{
					prefSize = std::max<size_t>(10, std::ceil(beamMaj*10.0));
					if(prefSize%2 != 0) ++prefSize;
				}
				++nIter;
			} while(!boxWasLargeEnough && nIter < 5);
		}
		else {
			fit2DGaussianCentred(image, width, height, beamEst, beamMaj, beamMin, beamPA);
		}
	}
	
private:
	const double* _image;
	size_t _width, _height, _scaleFactor;

	void fit2DGaussianCentredInBox(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA, size_t boxWidth, size_t boxHeight)
	{
		size_t startX = (width-boxWidth)/2;
		size_t startY = (height-boxHeight)/2;
		ao::uvector<double> smallImage(boxWidth*boxHeight);
		for(size_t y=startY; y!=(height+boxHeight)/2; ++y)
		{
			memcpy(&smallImage[(y-startY)*boxWidth], &image[y*width + startX], sizeof(double)*boxWidth);
		}
		
		fit2DGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamEst, beamMaj, beamMin, beamPA);
	}
	
	void fit2DGaussianCentred(const double* image, size_t width, size_t height, double beamEst, double& beamMaj, double& beamMin, double& beamPA)
	{
		_width = width;
		_height = height;
		_image = image;
		_scaleFactor = (width + height)/2;
		
		const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
		gsl_multifit_fdfsolver *solver = gsl_multifit_fdfsolver_alloc (T, _width*_height, 3);
		
		gsl_multifit_function_fdf fdf;
		fdf.f = &fitting_func;
		fdf.df = &fitting_deriv;
		fdf.fdf = &fitting_both;
		fdf.n = _width*_height;
		fdf.p = 3;
		fdf.params = this;
		
		// Using the FWHM formula for a Gaussian:
		const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
		double initialValsArray[3] = {
			beamEst/(_scaleFactor*double(sigmaToBeam)),
			beamEst/(_scaleFactor*double(sigmaToBeam)),
			0.0
		};
		gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 3);
		gsl_multifit_fdfsolver_set (solver, &fdf, &initialVals.vector);

		int status;
		size_t iter = 0;
		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (solver);
			
			if(status)
				break;
			
			status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);
			
		} while (status == GSL_CONTINUE && iter < 500);
		
		double
			sx = gsl_vector_get (solver->x, 0),
			sy = gsl_vector_get (solver->x, 1),
			beta = gsl_vector_get (solver->x, 2),
			betaFact = 1.0 - beta*beta;
		double cov[4];
		cov[0] = sx*sx / betaFact;
		cov[1] = beta * sx*sy / betaFact;
		cov[2] = cov[1];
		cov[3] = sy*sy / betaFact;
		gsl_multifit_fdfsolver_free(solver);
		
		double e1, e2, vec1[2], vec2[2];
		Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
		beamMaj = sqrt(std::fabs(e1)) * sigmaToBeam * _scaleFactor;
		beamMin = sqrt(std::fabs(e2)) * sigmaToBeam * _scaleFactor;
		beamPA = atan2(vec1[0], vec1[1]);
	}
	
	static int fitting_func(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double sx = gsl_vector_get(xvec, 0);
		double sy = gsl_vector_get(xvec, 1);
		double beta = gsl_vector_get(xvec, 2);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0/fitter._scaleFactor;
		
		size_t dataIndex = 0;
		double errSum = 0.0;
		for(size_t yi=0; yi!=height; ++yi)
		{
			double y = (yi - yMid) * scale;
			for(size_t xi=0; xi!=width; ++xi)
			{
				double x = (xi - xMid) * scale;
				double e = err(fitter._image[dataIndex], x, y, sx, sy, beta);
				errSum += e*e;
				gsl_vector_set(f, dataIndex, e);
				++dataIndex;
			}
		}
		//std::cout << "sx=" << sx << ", sy=" << sy << ", beta=" << beta << ", err=" << errSum << '\n';
		return GSL_SUCCESS;
	}
	
	static int fitting_deriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		GaussianFitter& fitter=*static_cast<GaussianFitter*>(data);
		double sx = gsl_vector_get(xvec, 0);
		double sy = gsl_vector_get(xvec, 1);
		double beta = gsl_vector_get(xvec, 2);
		const size_t width = fitter._width, height = fitter._height;
		int xMid = width/2, yMid = height/2;
		double scale = 1.0 / fitter._scaleFactor;
		
		size_t dataIndex = 0;
		for(size_t yi=0; yi!=height; ++yi)
		{
			double y = (yi - yMid)*scale;
			for(size_t xi=0; xi!=width; ++xi)
			{
				double x = (xi - xMid)*scale;
				double expTerm = exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
				double dsx = (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * expTerm;
				double dsy = (beta*x*y/(sy*sy*sx)+y*y/(sy*sy*sy)) * expTerm;
				double dbeta = -x*y/(sx*sy) * expTerm;
				gsl_matrix_set(J, dataIndex, 0, dsx);
				gsl_matrix_set(J, dataIndex, 1, dsy);
				gsl_matrix_set(J, dataIndex, 2, dbeta);
				++dataIndex;
			}
		}
		return GSL_SUCCESS;
	}
	
	static int fitting_both(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func(x, data, f);
		fitting_deriv(x, data, J);
		return GSL_SUCCESS;
	}
	
	static double err(double val, double x, double y, double sx, double sy, double beta)
	{
		return exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy)) - val;
	}
	static double dfdsx(double x, double y, double sx, double sy, double beta)
	{
		// f[x,y,sx,sy,beta]:=exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
		// diff(f[x,y,sx,sy,beta],sx);
		return (beta*x*y/(sx*sx*sy)+x*x/(sx*sx*sx)) * exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
	static double dfdbeta(double x, double y, double sx, double sy, double beta)
	{
		return x*y/(sx*sy) * exp(-x*x/(2.0*sx*sx) - beta*x*y/(sx*sy) - y*y/(2.0*sy*sy));
	}
};

#endif
