/**
 * @file amplici_umi.h
 * @author Xiyu Peng
 *
 * Cluster amplicon sequences with Unique Molecular Identifier
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#ifndef __AMPLICIUMI_H__
#define __AMPLICIUMI_H__

#include "options.h"
#include "initialize.h"
#include "amplici.h"

typedef struct {
	double omega;
	double pho;
	double *xi;
	unsigned int *idx;
	int32_t signs;   // Is it long enough ?
	int n_xi;
} mstep_data ;

typedef double (*msp_lcstr_func)(double lambda, void *fdata);
typedef double (*msp_lcstr_derv)(double lambda, void *fdata);

double ini_eta_gamma(options *opt, initializer *ini, double *gamma, double *eta, size_t sample_size);
int EM_algorithm(options *opt, data *dat, model *mod, initializer *ini, run_info *ri);
int log_vector(double *x, unsigned int len);
int normalize(size_t n, unsigned int K, double * Emat);
double mstep_newton(double (*fx)(double x, void *data), 
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata);
double mstep_pen1_lagrange_cstr_derv(double lambda, void *fdata);
double mstep_pen1_lagrange_cstr_func(double lambda, void *fdata);


#endif
