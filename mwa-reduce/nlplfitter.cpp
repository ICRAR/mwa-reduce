#include "nlplfitter.h"

#include <stdexcept>
#include <cmath>
#include <limits>

#include <iostream>

#ifdef HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#endif

class NLPLFitterData
{
public:
	typedef std::vector<std::pair<double, double>> PointVec;
	PointVec points;
#ifdef HAVE_GSL
	gsl_multifit_fdfsolver *solver;
	
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
	
	static int fitting_2nd_order(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
		double exponent = gsl_vector_get(xvec, 0);
		double facB = gsl_vector_get(xvec, 1);
		double facC = gsl_vector_get(xvec, 2);
		
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			double
				x = fitterData.points[i].first,
				y = fitterData.points[i].second;
			
			// f(x) = (bx + cx^2)^a
			gsl_vector_set(f, i, pow(facB*x + facC*x*x, exponent) - y);
		}
			
		return GSL_SUCCESS;
	}
	
	static int fitting_2nd_order_deriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
		double a = gsl_vector_get(xvec, 0);
		double b = gsl_vector_get(xvec, 1);
		double c = gsl_vector_get(xvec, 2);
	
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			double
				x = fitterData.points[i].first;
			
			// f(x)    = (bx + cx^2)^a
			// f(x)/da = ln(bx + cx^2) (bx + cx^2)^a
			// f(x)/db =    ax (bx + cx^2)^(a-1)
			// f(x)/dc =  ax^2 (bx + cx^2)^(a-1)
			double innerTerm = b*x + c*x*x;
			double toTheE = pow(innerTerm, a);
			double dfdexp = log(innerTerm) * toTheE;
			double toTheEM1 = toTheE/innerTerm;
			double dfdfacB = a*x*toTheEM1;
			double dfdfacC = a*x*x*toTheEM1;
				
			gsl_matrix_set(J, i, 0, dfdexp);
			gsl_matrix_set(J, i, 1, dfdfacB);
			gsl_matrix_set(J, i, 2, dfdfacC);
		}
			
		return GSL_SUCCESS;
	}

	static int fitting_2nd_order_both(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_2nd_order(x, data, f);
		fitting_2nd_order_deriv(x, data, J);
		return GSL_SUCCESS;
	}
#endif
};

#ifdef HAVE_GSL
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
	
	gsl_multifit_fdfsolver_free(_data->solver);
}

void NonLinearPowerLawFitter::Fit(double& a, double& b, double& c)
{
	Fit(a, b);
	b = pow(b, 1.0/a);
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	_data->solver = gsl_multifit_fdfsolver_alloc (T, _data->points.size(), 3);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &NLPLFitterData::fitting_2nd_order;
	fdf.df = &NLPLFitterData::fitting_2nd_order_deriv;
	fdf.fdf = &NLPLFitterData::fitting_2nd_order_both;
	fdf.n = _data->points.size();
	fdf.p = 3;
	fdf.params = &*_data;
	
	double initialValsArray[3] = { a, b, c };
	gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 3);
	gsl_multifit_fdfsolver_set(_data->solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (_data->solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(_data->solver->dx, _data->solver->x, 1e-7, 1e-7);
		
  } while (status == GSL_CONTINUE && iter < 500);
	
	a = gsl_vector_get (_data->solver->x, 0);
	b = gsl_vector_get (_data->solver->x, 1);
	c = gsl_vector_get (_data->solver->x, 2);
	
	gsl_multifit_fdfsolver_free(_data->solver);
}

#else
#warning "No GSL found: can not do non-linear power law fitting!"

void NonLinearPowerLawFitter::Fit(double& exponent, double& factor)
{
	throw std::runtime_error("Non-linear power law fitter was invoked, but GSL was not found during compilation, and is required for this");
}

#endif

NonLinearPowerLawFitter::NonLinearPowerLawFitter() :
	_data(new NLPLFitterData())
{
}

NonLinearPowerLawFitter::~NonLinearPowerLawFitter()
{
}

void NonLinearPowerLawFitter::AddDataPoint(double x, double y)
{
	_data->points.push_back(std::make_pair(x, y));
}

void NonLinearPowerLawFitter::FastFit(double& exponent, double& factor)
{
	double sumxy = 0.0, sumx = 0.0, sumy = 0.0, sumxx = 0.0;
	size_t n = 0;
	bool requireNonLinear = false;
	
	for(NLPLFitterData::PointVec::const_iterator i=_data->points.begin(); i!=_data->points.end(); ++i)
	{
		double x = i->first, y = i->second;
		if(y <= 0)
		{
			requireNonLinear = true;
			break;
		}
		if(x > 0 && y > 0)
		{
			long double
				logx = std::log(x),
				logy = std::log(y);
			sumxy += logx * logy;
			sumx += logx;
			sumy += logy;
			sumxx += logx * logx;
			++n;
		}
	}
	if(requireNonLinear)
	{
		exponent = 0.0;
		factor = 1.0;
		Fit(exponent, factor);
	}
	else {
		if(n == 0)
		{
			exponent = std::numeric_limits<double>::quiet_NaN();
			factor = std::numeric_limits<double>::quiet_NaN();
		}
		else {
			exponent = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
			factor = std::exp((sumy - exponent * sumx) / n);
		}
	}
}
