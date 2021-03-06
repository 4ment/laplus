#include "laplus.h"

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


double laplus_minimize(gsl_function* F, double guess, double lower, double upper, size_t max_iter){
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

void laplus_gamma_fit(double(*logP)(void* data, size_t, double),
						void(*dlogPs)(void* data, size_t, double*, double*),
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
		double dlogP, d2logP;
		dlogPs(data, i, &dlogP, &d2logP);
		
		double scale = -1.0/(map*d2logP);
		double shape = map/scale + 1.0;
		
		// Very small branch -> exponential shape
		if (map < 1.e-6 || d2logP >= 0) {
			shape = 1.0;
			scale = 1.0/fabs(dlogP);
			data_brent.value = scale;
			
			log_spaced_spaced_vector(x, map, 0.5, N);
			for (size_t j = 1; j < N; j++) {
				y[j] = logP(data, i, x[j]);
			}
			y[0] = logP(data, i, map); // leaves the tree with its original branch
			
			double maxY = y[0];
			for (size_t j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			
			gsl_function F;
			F.function = &func_gamma_fixed_shape;
			F.params = &data_brent;
			
			double guess = 1.0 - 0.001;
			if(func_gamma_fixed_shape(0.001, &data_brent) > func_gamma_fixed_shape(guess, &data_brent)
			   && func_gamma_fixed_shape(guess, &data_brent) < func_gamma_fixed_shape(1, &data_brent)){
				shape = laplus_minimize(&F, guess, 0.001, 1, 1000);
			}
		}
		// Small branch with a maximum and spurious large variance
		else if(shape*(scale*scale) > 0.1 && map < 0.0001){
			shape = 1.0;
			scale = 1.0/fabs(dlogP);
			data_brent.value = map;
			
			log_spaced_spaced_vector(x, map, 0.5, N);
			for (size_t j = 1; j < N; j++) {
				y[j] = logP(data, i, x[j]);
			}
			y[0] = logP(data, i, map); // leaves the tree with its original branch
			
			double maxY = y[0];
			for (size_t j = 0; j < N; j++) {
				y[j] -= maxY;
			}
			
			gsl_function F;
			F.function = &func_gamma_fixed_mode;
			F.params = &data_brent;
			
			double guess = 1.0 + 0.001;
			if(func_gamma_fixed_mode(1, &data_brent) > func_gamma_fixed_mode(guess, &data_brent)
			   && func_gamma_fixed_mode(guess, &data_brent) < func_gamma_fixed_mode(100, &data_brent)){
				shape = laplus_minimize(&F, guess, 1, 100, 1000);
				scale = map/(shape - 1);
			}
		}
	
		parameters[i*2] = shape;
		parameters[i*2+1] = scale;
	}
	free(x);
	free(y);
	free(yy);
}

double laplus_gamma_with_parameters(double map, const double* mapes, size_t count, const double* parameters){
	double  logP = map;
	for (size_t i = 0; i < count; ++i) {
		logP -= log(gsl_ran_gamma_pdf(mapes[i], parameters[i*2], parameters[i*2+1]));
	}
	return logP;
}

double laplus_gamma(double(*logP)(void* data, size_t, double),
					  void(*dlogPs)(void* data, size_t, double*, double*),
					  void* data,
					  size_t count,
					  double map, const double* mapes){
	
	double* parameters = calloc(count*2, sizeof(double));
	laplus_gamma_fit(logP, dlogPs, data, mapes, count, parameters);
	double logMarginal = laplus_gamma_with_parameters(map, mapes, count, parameters);
	free(parameters);
	return logMarginal;
}
