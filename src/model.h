/**
 * @file model.h
 * @author Karin S. Dorman
 *
 * Header file for model struct.
 */

#ifndef __H_MODEL__
#define __H_MODEL__

#include "data.h"
#include "options.h"

typedef struct _model model;	/* [KSD: move to model.h] */

/**
 * Parameter sets. [KSD: move to model.h]
 */
enum {
	PARAM_HAPLOTYPE = 1,
	PARAM_DELTA = 2,
	PARAM_LAMBDA = 4,
	PARAM_GAMMA = 8,
	PARAM_PI = 16,
	PARAM_BG_PI = 32,
};

/**
 * Model for bins
 * */
enum {
	NO_BINS = 0,
	BINS_LAMBDA0 = 1,
  	BINS_LAMBDA1 = 2,
	BINS_QUALITY = 4,
};

enum{
	EQUAL_LENGTH = 2,
	EXPECTATION_MODEL = 6,  /* for lambda0 only */
};
 

struct _model {

	unsigned char n_quality;/*<! no. of distinct quality scores after compression */

	double *eik;		/*<! expectations or transition prob (used in UMI model) */
	unsigned int n_param;	/*<! number of parameters */
	unsigned int K; /*<! number of clusters */
	double *pi;		/*<! kx1 mixing proportions */
	unsigned char *haplotypes;	/*<! haplotypes */
	double ll;		/*<! current log likelihood */
	double best_ll;		/*<! best log likelihood so far */
	double JC_ll;      /*<! log likelihood of haplotypes of JC69 model */
	
	/* error profile */
	double *error_profile;		/*<! error profile */
	int err_encoding;      /*<! nucleotide encodings of error profile */
	double adj_trunpois; /*<! adjust value for truncated poisson */

	/* model comparison */
	double aic;
	double bic;

	/* for the JC69 model */
	double * distance; /*<! Kx1 distances from haplotypes to the ancestor */
	unsigned char * est_ancestor; /*<! The estimated ancestor */
	double *JC_ll_K; /*<! Kx1 log likelihood under JC69 model */

	/* for UMI model */
	double *eik_umi;   /*<! transition prob of each observed UMIs */
	double *gamma;    /*<! K_umi x K transiton prob of UMIs to haplotypes (log) */
	double *eta;      /*<! K_umi proportions of UMIs (log) */
	unsigned int *E2_sparse_hap_id; /*<! store Hap_id (topN x N) */
	unsigned int *E2_sparse_umi_id; /*<! store UMI_id (topN x N) */
	double *E2_sparse_value; /*<! store posterial prob of reads (topN x N) */
	double ll_UMI;    /*<! current log likelihood for UMI model */
	double pll_UMI;    /*<! current log likelihood for UMI model (previous copy) */
	double penalty_ll;   /*<! penalty for transition prob */


}; /* model */

double translate_error_STD_to_XY(double *error_profile, unsigned char n_quality,
	unsigned char hap_nuc, unsigned char obser_nuc, unsigned char qual);

int make_model(model **mod, data *dat, options *opt);
int realloc_model(model *mod, data *dat, options *opt);
void free_model(model *mod);

#endif
