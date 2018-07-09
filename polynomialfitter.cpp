#include "polynomialfitter.h"

#include <gsl/gsl_multifit.h>

void PolynomialFitter::Fit(ao::uvector<double>& terms, size_t nTerms)
{
	size_t n = _dataPoints.size();
	terms.assign(nTerms, 0.0);
	
	if(nTerms > n)
	{
		nTerms = n;
	}
	
	gsl_multifit_linear_workspace*
		workspace = gsl_multifit_linear_alloc(n, nTerms);
		
	gsl_matrix* xData = gsl_matrix_alloc(n, nTerms);
	gsl_matrix* cov = gsl_matrix_alloc(nTerms, nTerms);
	gsl_vector* wData = gsl_vector_alloc(n);
	gsl_vector* yData = gsl_vector_alloc(n);
	gsl_vector* resultTerms = gsl_vector_alloc(nTerms);
	double chisq;
	
	for(size_t i=0; i!=n; ++i)
	{
		double x = _dataPoints[i][0];
		double y = _dataPoints[i][1];
		double w = _dataPoints[i][2];
		
		double f = 1.0;
		gsl_matrix_set(xData, i, 0, f);
		for(size_t j=1; j!=nTerms; ++j)
		{
			f *= x; // f = x^j
			gsl_matrix_set(xData, i, j, f);
		}
		gsl_vector_set(yData, i, y);
		gsl_vector_set(wData, i, w);
	}
	
	int result = gsl_multifit_wlinear(xData, wData, yData, resultTerms, cov, &chisq, workspace);
	
	for(size_t j=0; j!=nTerms; ++j)
	{
		terms[j] = gsl_vector_get(resultTerms, j);
	}

	gsl_vector_free(resultTerms);
	gsl_vector_free(yData);
	gsl_vector_free(wData);
	
	gsl_matrix_free(cov);
	gsl_matrix_free(xData);
	
	gsl_multifit_linear_free(workspace);
	
	if(result)
	{
		throw std::runtime_error("Polynomial fit failed");
	}
}
