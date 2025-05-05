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
	options *opt;			/*<! pointer to options object */

	unsigned char n_quality;	/*<! no. post-compression quality scores */

	double *eik;			/*<! expectated no. seq i from hap k */
	unsigned int K;			/*<! number of clusters */
	double *pi;			/*<! Kx1 mixing proportions */
	unsigned char *haplotypes;	/*<! haplotypes */
	double ll;			/*<! current log likelihood */
	double best_ll;			/*<! best log likelihood so far */
	double JC_ll;      		/*<! log likelihood of JC69 model */
	
	/* error model */
	double *error_profile;		/*<! error profile */
	double *precomputed_dindel;	/*<! precomputed indel distn (fixed version) */
	double p_indel;			/*<! prob indel: from options::indel_error */
	int err_encoding;		/*<! nucleotide order in error profile */
	unsigned int rd_length;		/*<! read length */
		/* (OG AmpliCI on variable-length reads effectively set to max) */
	double adj_trunpois;		/*<! Pr(X<=t) for truncated poisson */

	/* model comparison */
	double aic;
	double bic;
	unsigned int n_param;		/*<! number of parameters */

	/* for the JC69 model */
	double *distance;		/*<! Kx1 haplotypes to ancestor distances */
	unsigned char *est_ancestor;	/*<! estimated ancestor */
	double *JC_ll_K;		/*<! Kx1 log likelihood under JC69 model */

	/* for UMI model */
	double *eik_umi;		/*<! transition prob of observed UMIs */
	double *gamma;			/*<! K_umi x K log UMIs to hap trans. prob */
	double *eta;			/*<! K_umi proportions of UMIs (log) */
	unsigned int *E2_sparse_hap_id;	/*<! store Hap_id (topN x N) */
	unsigned int *E2_sparse_umi_id;	/*<! store UMI_id (topN x N) */
	double *E2_sparse_value;	/*<! store posterior prob of reads (topN x N) */
	double ll_UMI;			/*<! current log likelihood */
	double pll_UMI;			/*<! previous log likelihood */
	double penalty_ll;		/*<! penalty for transition prob */


}; /* model */

double translate_error_STD_to_XY(double *error_profile, unsigned char n_quality, unsigned char hap_nuc, unsigned char obser_nuc, unsigned char qual);
double dindel(model *mod, unsigned int n_indel, unsigned int n_opp, int logged);

int make_model(model **mod, data *dat, options *opt);
int realloc_model(model *mod, data *dat, options *opt);
void free_model(model *mod);

#endif
