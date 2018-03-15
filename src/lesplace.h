//
//  lesplace.h
//  lesplace
//
//  Created by Mathieu Fourment on 15/3/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#ifndef lesplace_h
#define lesplace_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void lesplace_gamma_fit(double(*likelihood)(void* data, size_t, double), // returns the log joint
					double(*dlikelihood)(void* data, size_t), // returns the first derivative of the log joint wrt the indexed parameter
					double(*d2likelihood)(void* data, size_t), // returns the second derivative of the log joint wrt the indexed parameter
					void(*reset)(void* data, size_t, double), // reset the value of the indexed parameter
					void* data, // data for callback functions
					const double* mapes, size_t count, // maximum a posteriori estimates
					double* parameters); // parameters of the distribution
	
	// map: maximum a posteriori
	// mapes: maximum a posteriori estimates
	// parameters: parameters of the gamma distribution. Should be of size 2*count
double lesplace_gamma(double map, const double* mapes, size_t count, const double* parameters);
	
#ifdef __cplusplus
} // extern "C"
#endif

#endif /* lesplace_h */
