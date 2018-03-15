//
//  lesplace.c
//  lesplace
//
//  Created by Mathieu Fourment on 15/3/18.
//  Copyright © 2018 Mathieu Fourment. All rights reserved.
//

#include "lesplace.h"

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_min.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>


void log_spaced_spaced_vector( double *v, double lower, double upper, size_t N ){
	double logMin = log(lower);
	double logMax = log(upper);
	double delta = (logMax - logMin) / (N-1);
	
	double accDelta = 0;
	for (size_t i = 0; i < N; ++i) {
		v[i] = exp(logMin + accDelta);
		accDelta += delta;
	}
}

struct marginal_data_t{
	double* x;
	double* y;
	double*yy;
	size_t N;
	double value;
};

// shape is fixed and scale is optimized
double func_gamma_fixed_shape(double x, void* params){
	struct marginal_data_t* d = (struct marginal_data_t*)params;
	double scale = d->value;
	double maxY = -INFINITY;
	
	for (size_t i = 0; i < d->N; i++) {
		d->yy[i] = gsl_ran_gamma_pdf(d->x[i], x, scale);
		if(d->yy[i] > maxY){
			maxY = d->yy[i];
		}
	}
	double sum = 0;
	for (size_t i = 0; i < d->N; i++) {
		sum += pow((d->yy[i]-maxY) - d->y[i], 2);
	}
	return sum;
}

// scale is optimized and shape is constrained
double func_gamma_fixed_mode(double x, void* params){
	struct marginal_data_t* d = (struct marginal_data_t*)params;
	double shape = x;
	// mode = (shape-1)scale for shape >= 1
	double scale = d->value/(shape - 1.0);
	double maxY = -INFINITY;
	
	for (size_t i = 0; i < d->N; i++) {
		d->yy[i] = gsl_ran_gamma_pdf(d->x[i], shape, scale);
		if(d->yy[i] > maxY){
			maxY = d->yy[i];
		}
	}
	double sum = 0;
	for (size_t i = 0; i < d->N; i++) {
		sum += pow((d->yy[i]-maxY) - d->y[i], 2);
	}
	return sum;
}


double lesplace_minimize(gsl_function* F, double guess, double lower, double upper, size_t max_iter){
	const gsl_min_fminimizer_type* T = gsl_min_fminimizer_brent;
	gsl_min_fminimizer* s = gsl_min_fminimizer_alloc (T);
	gsl_min_fminimizer_set (s, F, guess, lower, upper);
	
	size_t iter = 0;
	double a, b, m;
	int status;
	
	do{
		iter++;
		status = gsl_min_fminimizer_iterate (s);
		
		m = gsl_min_fminimizer_x_minimum (s);
		a = gsl_min_fminimizer_x_lower (s);
		b = gsl_min_fminimizer_x_upper (s);
		
		status = gsl_min_test_interval (a, b, 0.001, 0.0);
		
		if (status == GSL_SUCCESS){
			gsl_min_fminimizer_free (s);
			return m;
		}
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	
	gsl_min_fminimizer_free (s);
	return INFINITY;
}

void lesplace_gamma_fit(double(*likelihood)(void* data, size_t, double),
						double(*dlikelihood)(void* data, size_t),
						double(*d2likelihood)(void* data, size_t),
						void(*reset)(void* data, size_t, double),
						void* data,
						const double* mapes, size_t count,
						double* parameters){
	
	size_t N = 10;
	double* x = calloc(N, sizeof(double));
	double* y = calloc(N, sizeof(double));
	double* yy = calloc(N, sizeof(double));
	struct marginal_data_t data_brent = {x, y, yy, N, 0};
	
	for (size_t i = 0; i < count; ++i) {
		const double map = mapes[i];
		double d2logP = d2likelihood(data, i);
		
		double scale = -1.0/(map*d2logP);
		double shape = map/scale + 1.0;
		
		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			double dlogP = dlikelihood(data, i);
			scale = 1.0/fabs(dlogP);
			data_brent.value = scale;
			
			log_spaced_spaced_vector(x, map, 0.5, N);
			for (size_t j = 0; j < N; j++) {
				y[j] = likelihood(data, i, x[j]);
			}
			double maxY = y[0];
			for (size_t j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			
			gsl_function F;
			F.function = &func_gamma_fixed_shape;
			F.params = &data_brent;
			
			shape = lesplace_minimize(&F, 1, 0, 1, 1000);
			
			reset(data, i, map);
		}
		// Small branch with a maximum and spurious large variance
		else if(shape/(scale*scale) > 0.1 && map < 0.0001){
			data_brent.value = map;
			
			log_spaced_spaced_vector(x, map, 0.5, N);
			for (size_t j = 0; j < N; j++) {
				y[j] = likelihood(data, i, x[j]);
			}
			double maxY = y[0];
			for (size_t j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			
			gsl_function F;
			F.function = &func_gamma_fixed_mode;
			F.params = &data_brent;
			
			shape = lesplace_minimize(&F, 2, 1, 100, 1000);
			scale = map/(shape - 1);
			
			reset(data, i, map);
		}
	
		parameters[i*2] = shape;
		parameters[i*2+1] = scale;
	}
	free(x);
	free(y);
	free(yy);
}

double lesplace_gamma(double map, const double* mapes, size_t count, const double* parameters){
	double  logP = map;
	for (size_t i = 0; i < count; ++i) {
		logP -= log(gsl_ran_gamma_pdf(mapes[i], parameters[i*2], parameters[i*2+1]));
	}
	return logP;
}

