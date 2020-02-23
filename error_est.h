/**
 * @file error_est.c
 * @author Xiyu Peng
 * 
 * Header file for error_est.h
 */

#ifndef __H_ERROR_EST__
#define __H_ERROR_EST__

#include "amplici.h"
#include "ampliclust.h"
#include "statistics.h"
#include "lmath.h"
#include "io.h"
#include "error.h"
#include "loess.h"

typedef struct {
    double beta0;
    double beta1;
} linear;

int error_profile_generator(options *opt, data *dat, model *mod, initializer *ini, run_info *ri);
int error_predict(unsigned int *err_cnt, double *error_profile, unsigned int n_quality,unsigned int self_lines[4]);
void fprint_error_profile(FILE *fp, double *mat, unsigned int n, unsigned int l);

/* functions of regression or to call regression */
int loess_est(loess *lo,double *y,double *x, unsigned int* weights, unsigned int n);

/* for loess regression */
void loess_setup(double *x, double *y, int n, int p, loess *lo);
void loess_summary(loess *lo);
void loess_fit(loess *lo);
void loess_free_mem(loess *lo);

void predict(double *eval, int m, loess *lo, prediction *pre, int se);
void pred_free_mem(prediction *pre); 

double cos_dist(unsigned int *mat1, unsigned int *mat2, unsigned int nrow, unsigned int ncol);

#endif
