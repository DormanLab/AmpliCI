/**
 * @file options.h
 * @author Karin S. Dorman
 *
 * Header file for options struct.
 */

#ifndef __H_OPTIONS__
#define __H_OPTIONS__

#ifdef USE_CURSES
#include <curses.h>
#endif

#include "fastq.h"

enum {
	ALGORITHM_AMPLICLUST,   
	ALGORITHM_AMPLICI
};


typedef struct _options options;
/**
 * Run options.
 */
struct _options {
	
	/* ampliCI */
	int run_amplici;  /* if we only run ampliCI */
	double low_bound;  /*<! minimum haplotype abundance */
	/* input */
	char const *offset_file;  /*<! name of offset file */
	char const *fastq_file;  /*<! name of fastq input file */
	char const *initialization_file;   /*<! name of initialization file */
	/* output */
	char const *outfile;  /*<! ... */

	/* model */
	int convergence_amplici;  /*<! convergence or not when updating abundance */
	int check_false_positive;  /*<! if we check false positive */
	int use_aic;  /*<! use aic to do model selection */
	int JC69_model;  /*<! use approx. JC69 hierarchical model */
	int nw_align;  /*<! when and how we perform nw alignment */
	int indel_model;  /*<! build indel model to calculate trans prob*/

	/* error profile estimation */
	int error_estimation;  /*<! run error estimate function */
	int use_error_profile;  /*<! if we choose to use an input error profile */
	int err_encoding;  /*<! how nucleotides encoded in error profile */
	char const *error_profile_name;  /*<! name of file name of error profile */
	double min_cosdist;  /*<! minimum cosine distance (log version) */

	/* error */
	double insertion_error;  /*<! insertion error rate */
	double deletion_error;  /*<! deletion error rate */
	double indel_error;  /*<! indel error rate per site per read */
	
	/* K number of clusters */
	int estimate_K;  /*<! whether algorithm should estimate K */
	unsigned int K;  /*<! number of clusters */
	unsigned int K_max;  /*<! only select first K_max haplotypes */
	int K_fix_err;  /*<! fix the maximum K when estimate error profiles */
	unsigned int K_space;  /*<! record the size of memory */

	/* running criterion */
	unsigned int n_iter_amplici;  /*<! max. number of AmpliCI iterations */
	double ll_cutoff;  /*<! set a cutoff of the likelihood filter */
	double epsilon_aln;  /*<! trans prob when nw align is out of band  */
	double epsilon;  /*<! change in relative abundance */
	double alpha;  /*<! FWER */
	double p_threshold;  /*<! threshold on p-value */
	unsigned int most_abundant;  /*<! report these most abundant only */
	
	/* alignment */
	int score[NUM_NUCLEOTIDES][NUM_NUCLEOTIDES];  /* score matrix for nw alignment */
	int gap_p;  /* penalty for gaps */
	int band;  /* bandwidth for banded nw alignment */

	int info;  /*<! level of information to output */
	int use_curses;
#ifdef USE_CURSES
	WINDOW *wp;
#else
	FILE *wp;
#endif
	FILE *active_fp;
}; /* options */

int make_options(options **opt);
int parse_options(options *opt, int argc, const char **argv);
void free_options(options *opt);

#endif
