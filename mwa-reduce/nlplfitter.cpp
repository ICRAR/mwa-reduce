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
	typedef ao::uvector<std::pair<double, double>> PointVec;
	PointVec points;
	size_t nTerms;
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
	
	static int fitting_multi_order(const gsl_vector *xvec, void *data, gsl_vector *f)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
		
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			const double
				x = fitterData.points[i].first,
				y = fitterData.points[i].second;
			
			const double lg = log10(x);
			
			// Horner's method
			double fity = 0.0;
			for(size_t k=0; k!=fitterData.nTerms; ++k)
			{
				size_t j = fitterData.nTerms-k-1;
				const double a_j = gsl_vector_get(xvec, j);
				fity = a_j + fity * lg;
			}
			//std::cout << x << ':' << fity << " / " << exp10(fity) << "\n";
			gsl_vector_set(f, i, exp10(fity) - y);
		}
			
		return GSL_SUCCESS;
	}
	
	static int fitting_multi_order_deriv(const gsl_vector *xvec, void *data, gsl_matrix *J)
	{
		const NLPLFitterData &fitterData = *reinterpret_cast<NLPLFitterData*>(data);
	
		for(size_t i=0; i!=fitterData.points.size(); ++i)
		{
			double
				x = fitterData.points[i].first;
			
			const double lg = log10(x);
			
			// Horner's method
			double fity = 0.0;
			for(size_t k=0; k!=fitterData.nTerms; ++k)
			{
				size_t j = fitterData.nTerms-k-1;
				const double a_j = gsl_vector_get(xvec, j);
				fity = a_j + fity * lg;
			}
			fity = exp10(fity);
			// dY/da_i = e^[ a_0...a_i-1,a_i+1...a_n] * (e^[a_i {log x}^i]) {log x}^i
			gsl_matrix_set(J, i, 0, fity);
			
			double lgPower = lg;
			for(size_t j=1; j!=fitterData.nTerms; ++j)
			{
				//const double a_j = gsl_vector_get(xvec, j);
				gsl_matrix_set(J, i, j, fity*lgPower);
				lgPower *= lg;
			}
		}
			
		return GSL_SUCCESS;
	}

	static int fitting_multi_order_both(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
	{
		fitting_multi_order(x, data, f);
		fitting_multi_order_deriv(x, data, J);
		return GSL_SUCCESS;
	}
	
#endif
};

#ifdef HAVE_GSL
void NonLinearPowerLawFitter::Fit(double& exponent, double& factor)
{
	if(_data->points.size() >= 2)
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
	else {
		exponent = 0.0;
		factor = 0.0;
		for(size_t i=0; i!=_data->points.size(); ++i)
			factor += _data->points[i].second;
		factor /= double(_data->points.size());
	}
}

void NonLinearPowerLawFitter::Fit(double& a, double& b, double& c)
{
	Fit(a, b);
	b = pow(b, 1.0/a);
	
	if(_data->points.size() >= 3)
	{
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
}

void NonLinearPowerLawFitter::fit_implementation(ao::uvector<double>& terms, size_t nTerms)
{
	_data->nTerms = nTerms;
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
	_data->solver = gsl_multifit_fdfsolver_alloc (T, _data->points.size(), nTerms);
	
	gsl_multifit_function_fdf fdf;
	fdf.f = &NLPLFitterData::fitting_multi_order;
	fdf.df = &NLPLFitterData::fitting_multi_order_deriv;
	fdf.fdf = &NLPLFitterData::fitting_multi_order_both;
	fdf.n = _data->points.size();
	fdf.p = nTerms;
	fdf.params = &*_data;
	
	gsl_vector_view initialVals = gsl_vector_view_array(terms.data(), nTerms);
	gsl_multifit_fdfsolver_set(_data->solver, &fdf, &initialVals.vector);

	int status;
	size_t iter = 0;
	do {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (_data->solver);
		
		if(status)
			break;
		
		status = gsl_multifit_test_delta(_data->solver->dx, _data->solver->x, 1e-6, 1e-6);
		
  } while (status == GSL_CONTINUE && iter < 5000);
	
	if(status != GSL_SUCCESS)
	{
		std::cout << "Warning: not converged! (niter=" << iter << ", status=" << gsl_strerror(status) << ")\n";
		//for(size_t i=0; i!=nTerms; ++i)
		//	terms[i] = std::numeric_limits<double>::quiet_NaN();
	}
	for(size_t i=0; i!=nTerms; ++i)
		terms[i] = gsl_vector_get (_data->solver->x, i);
	
	gsl_multifit_fdfsolver_free(_data->solver);
}

#else
#warning "No GSL found: can not do non-linear power law fitting!"

void NonLinearPowerLawFitter::Fit(double& exponent, double& factor)
{
	throw std::runtime_error("Non-linear power law fitter was invoked, but GSL was not found during compilation, and is required for this");
}

void NonLinearPowerLawFitter::Fit(double& a, double& b, double& c)
{
	throw std::runtime_error("Non-linear power law fitter was invoked, but GSL was not found during compilation, and is required for this");
}

void NonLinearPowerLawFitter::fit_implementation(std::vector<double>& terms, size_t nTerms)
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

void NonLinearPowerLawFitter::Fit(ao::uvector<double>& terms, size_t nTerms)
{
	terms.assign(nTerms, 0.0);
	
	if(_data->points.size() < nTerms)
		nTerms = _data->points.size();
	
	if(nTerms == 0)
		return;
	
	double a = 1.0, b = 0.0;
	Fit(a, b);
	bool isNegative = b < 0.0;
	if(isNegative)
	{
		for(NLPLFitterData::PointVec::iterator i=_data->points.begin(); i!=_data->points.end(); ++i)
		{
			i->second = -i->second;
		}
		terms[0] = -log10(-b);
		a = -a;
	}
	else {
		terms[0] = log10(b); // - a*log(NLPLFact);
	}
	
	if(b != 0.0)
	{
		if(nTerms > 1) terms[1] = a;
		
		fit_implementation(terms, nTerms);
	}
	
	if(isNegative)
		terms[0] = -exp10(terms[0]);
	else
		terms[0] = exp10(terms[0]);
}

void NonLinearPowerLawFitter::FitStable(ao::uvector<double>& terms, size_t nTerms)
{
	terms.assign(nTerms, 0.0);
	if(nTerms == 0)
		return;
	
	double a = 1.0, b = 0.0;
	Fit(a, b);
	
	bool isNegative = b < 0.0;
	if(isNegative)
	{
		for(NLPLFitterData::PointVec::iterator i=_data->points.begin(); i!=_data->points.end(); ++i)
		{
			i->second = -i->second;
		}
		terms[0] = -log10(-b);
		a = -a;
	}
	else {
		terms[0] = log10(b); // - a*log(NLPLFact);
	}
		
	if(b != 0.0)
	{
		if(nTerms > 1) terms[1] = a;
		size_t nTermsEstimated = 2;
		while(nTermsEstimated < nTerms)
		{
			++nTermsEstimated;
			fit_implementation(terms, nTermsEstimated);
		}
	}
	
	if(isNegative)
		terms[0] = -exp10(terms[0]);
	else
		terms[0] = exp10(terms[0]);
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
			double d = (n * sumxx - sumx * sumx);
			if(d == 0.0)
				exponent = 0.0;
			else
				exponent = (n * sumxy - sumx * sumy) / d;
			factor = std::exp((sumy - exponent * sumx) / n);
		}
	}
}

double NonLinearPowerLawFitter::Evaluate(double x, const ao::uvector<double>& terms, double referenceFrequencyHz)
{
	if(terms.empty()) return 0.0;
	double y = 0.0;
	const double lg = log10(x/referenceFrequencyHz);
	
	for(size_t k=0; k!=terms.size()-1; ++k)
	{
		size_t j = terms.size()-k-1;
		y = y * lg + terms[j];
	}
	return exp10(y * lg) * terms[0];
}
