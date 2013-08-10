#include "nlplfitter.h"
#include <stdexcept>

#ifdef HAVE_GSL

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

class NLPLFitterData
{
public:
	typedef std::vector<std::pair<double, double>> PointVec;
	gsl_multifit_fdfsolver *solver;
	PointVec points;
	
	static int fitting_func(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
		double exponent = gsl_vector_get(xvec, 0);
		double factor = gsl_vector_get(xvec, 1);
		
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			double
				x = fitterData.points[i].first,
				y = fitterData.points[i].second;
			
			gsl_vector_set(f, i, factor * pow(x, exponent) - y);
		}
			
		return GSL_SUCCESS;
	}
	
	static int fitting_func_deriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
		double exponent = gsl_vector_get(xvec, 0);
		double factor = gsl_vector_get(xvec, 1);
	
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			double
				x = fitterData.points[i].first;
			
			double xToTheE = pow(x, exponent);
			double dfdexp = factor * log(x) * xToTheE;
			double dfdfac = xToTheE;
				
			gsl_matrix_set(J, i, 0, dfdexp);
			gsl_matrix_set(J, i, 1, dfdfac);
		}
			
		return GSL_SUCCESS;
	}

	static int fitting_func_both(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_func(x, data, f);
		fitting_func_deriv(x, data, J);
		return GSL_SUCCESS;
	}
};

NonLinearPowerLawFitter::NonLinearPowerLawFitter() :
	_data(new NLPLFitterData())
{
}

NonLinearPowerLawFitter::~NonLinearPowerLawFitter()
{
}

void NonLinearPowerLawFitter::Fit(double& exponent, double& factor)
{
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	_data->solver = gsl_multifit_fdfsolver_alloc (T, _data->points.size(), 2);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &NLPLFitterData::fitting_func;
	fdf.df = &NLPLFitterData::fitting_func_deriv;
	fdf.fdf = &NLPLFitterData::fitting_func_both;
	fdf.n = _data->points.size();
	fdf.p = 2;
	fdf.params = &*_data;
	
	double initialValsArray[2] = { exponent, factor };
	gsl_vector_view initialVals = gsl_vector_view_array (initialValsArray, 2);
	gsl_multifit_fdfsolver_set (_data->solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (_data->solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(_data->solver->dx, _data->solver->x, 1e-7, 1e-7);
		
  } while (status == GSL_CONTINUE && iter < 500);
	
	exponent = gsl_vector_get (_data->solver->x, 0);
	factor = gsl_vector_get (_data->solver->x, 1);
}

void NonLinearPowerLawFitter::AddDataPoint(double x, double y)
{
	_data->points.push_back(std::make_pair(x, y));
}

#else
#warning "No GSL found: can not do non-linear power law fitting!"

class NLPLFitterData
{
};

NonLinearPowerLawFitter::NonLinearPowerLawFitter()
{
	throw std::runtime_error("Non-linear power law fitter was invoked, but GSL was not found during compilation, and is required for this");
}

NonLinearPowerLawFitter::~NonLinearPowerLawFitter()
{
}

void NonLinearPowerLawFitter::AddDataPoint(double x, double y)
{
}

void NonLinearPowerLawFitter::Fit(double& exponent, double& factor)
{
}

#endif
