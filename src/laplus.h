#ifndef laplus_h
#define laplus_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void laplus_gamma_fit(double(*logP)(void* data, size_t, double), // returns the log joint
					void(*dlogPs)(void* data, size_t, double*, double*), // returns the first and second derivatives of the log joint wrt the indexed parameter
					void* data, // data for callback functions
					const double* mapes, size_t count, // maximum a posteriori estimates
					double* parameters); // parameters of the distribution
	
	// map: maximum a posteriori
	// mapes: maximum a posteriori estimates
	// parameters: parameters of the gamma distribution. Should be of size 2*count
double laplus_gamma_with_parameters(double map, const double* mapes, size_t count, const double* parameters);
	
double laplus_gamma(double(*logP)(void* data, size_t, double),
					  void(*dlogPs)(void* data, size_t, double*, double*),
					  void* data,
					  size_t count,
					  double map, const double* mapes);
	
#ifdef __cplusplus
} // extern "C"
#endif

#endif /* laplus_h */
