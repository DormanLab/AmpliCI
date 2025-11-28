/**
 * @file amplici.h
 * @author Xiyu Peng
 *
 * Header file for amplici.c
 */

#ifndef __H_AMPLICI__
#define __H_AMPLICI__

#include "model.h"
#include "data.h"
#include "options.h"
#include "initialize.h"
#include "ampliclust.h"

enum {
	NO_ALIGNMENT,  /*<! model with no gap */
	ALIGNMENT_UNIQ_SEQ,    /*<! nw alignment between haplotype and uniq seq */
	ALIGNMENT_HAPLOTYPES    /*<! nw alignment between haplotypes */
};

enum {
	IGNORE_INDELS,    /*< ! only consider the overlapping part */
	INDEL_PER_SITE1,		/*< ! use specific indel error rate */
	INDEL_PER_SITE2,    /*<! calculate prob for every site of the read */
	INDEL_PER_READ		/* <! calculate indel error per read */
};

#define MAX_N_EXACT_P 100

/* main function */
int ampliCI(options * opt, data * dat, model *mod, initializer *ini, run_info *ri);
int haplotype_selection(options * opt, data * dat, model *mod, initializer *ini, unsigned int K_max);
int reads_assignment(options * opt, data * dat, model *mod, initializer *ini, run_info *ri);

/* others */
double Simple_Estep(model *mod, size_t sample_size, double *e_trans, unsigned int K);
int likelihood_filter(unsigned int K, double ll_cutoff, double *eik, double *pi, double * e_trans, size_t sample_size, run_info *ri);
int trans_expectation(options *opt, model *mod, data *dat, initializer*ini, double *error_profile, double *trans_prob, int ends_free);
double trans_nw(options *opt, model *mod, data_t **align, size_t alen, unsigned int mismatch, unsigned int ngap, double *error_profile, int err_encoding, data_t *rqmat, data_t n_quality, unsigned int rlen, double *error_prob, int ends_free);

/* print */
void fprint_assignment(FILE *fp, unsigned int *v, size_t n, unsigned int max, int width, int newline);
void fprint_haplotypes_abun(FILE *fp, data_t **data, size_t n, unsigned int *len, double pthres, char const * const prefix, double *pvalue, double *abun, double *ee);
void fprint_haplotypes_size(FILE *fp, data_t **data, size_t n, unsigned int *len, double pthres, char const * const prefix, double *pvalue, unsigned int *size, double *ee);

/* nwalign */
int nwalign_matrix(options *opt, data *dat, initializer *ini, data_t ***nw_result, size_t *nw_alen, unsigned int size, unsigned int select);
void free_nw_result(data_t ***nw_result, unsigned int space);
//void realloc_nw_result(data_t ***nw_result, unsigned int space);

/* sequence errors */
double exp_errors(data_t *qual, unsigned int length, double *error_prob);

/* BIC */
int m_JC69(data_t *hap, data_t *anc, double *dist, unsigned int K, unsigned int len, int start);
double e_JC69(data_t *hap, data_t *anc, double *dist, unsigned int K, unsigned int len, int start);
int modified_ic(data_t *hap, data_t *est_anc, double *distance, double best_ll, unsigned int K, double *JC_ll, double *n_aic, double *n_bic, unsigned int n_param, unsigned int max_read_length, size_t sample_size, int start);

#endif
