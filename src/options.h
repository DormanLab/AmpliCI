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

enum {
	MLE,
	MPLE
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
	char const *initialization_file;  /*<! name of initialization file */
	char const *trans_matrix; /*<! output transition matrix file */
	/* output */
	char const *outfile_base;  /*<! basename of outfiles */
	char const *outfile_info;  /*<! name of informational outfile */
	char const *outfile_fasta;  /*<! name of fasta outfile */

	/* UMI information */
	unsigned int UMI_length;  /* length of UMIs (at the beginning of reads ) */
	char const *initialization_UMI;  /*<! name of initialization file for UMIs */
	unsigned int K_UMI;  /*<! num of UMI clusters */
	unsigned int topN;  /*<! top N possible combinations we considered in the model */
	int trans_penalty;   /*<! penalty on transiton prob */
	double rho;          /*<! for MPLE */
	double omega;       /*< ! for MPLE */
	double threshold_UMI;   /* ! minimal UMI abundance */
	double threshold_hap;   /*  !minimal hap abundance */
	unsigned int max_offset; /* ! maximal allowed offset */
	unsigned int umicollision; /* ! consider UMI collision or not*/

	/* model */
	int convergence_amplici;  /*<! convergence or not when updating abundance */
	int check_false_positive;  /*<! if we check false positive */
	int use_aic;  /*<! use aic to do model selection */
	int JC69_model;  /*<! use approx. JC69 hierarchical model */
	int nw_align;  /*<! when and how we perform nw alignment */
	int indel_model;  /*<! build indel model to calculate trans prob*/
	int ignor_nc; /*<! Number of nucleotides that will be ignored in JC69 hierarchical model*/

	/* error profile estimation */
	int error_estimation;  /*<! run error estimate function */
	char const *partition_file; /*<! given partition file */
	unsigned int seed_min_observed_abundance; /*<! min. abund. for seed */
	unsigned char exclude_low_abundance_seeds; /*<! exclude low abundance seeds */
	int use_error_profile;  /*<! if we choose to use an input error profile */
	int err_encoding;  /*<! how nucleotides encoded in error profile */
	char const *error_profile_name;  /*<! name of file name of error profile */
	double min_cosdist;  /*<! minimum cosine distance (log version) */

	/* error */
	double insertion_error;  /*<! insertion error rate */
	double deletion_error;  /*<! deletion error rate */
	double indel_error;  /*<! indel error rate per site per read */
	int indel_error_set;  /*<! has indel error been set? */
	
	/* K number of clusters */
	int estimate_K;  /*<! whether algorithm should estimate K */
	unsigned int K;  /*<! number of clusters */
	unsigned int K_max;  /*<! only select first K_max haplotypes */
	int K_fix_err;  /*<! fix the maximum K when estimate error profiles */
	unsigned int K_space;  /*<! record the size of memory */
	int filter_reads; /*<! filter reads in error estimation */


	/* running criterion */
	unsigned int n_iter_amplici;  /*<! max. number of AmpliCI iterations */
	double ll_cutoff;  /*<! set a cutoff of the likelihood filter */
	double epsilon_aln;  /*<! trans prob when nw align is out of band  */
	double epsilon;  /*<! change in relative abundance */
	double alpha;  /*<! FWER */
	double p_threshold;  /*<! threshold on p-value */
	int per_candidate;  /*<! threshold on p-value is per candidate */
	unsigned int most_abundant;  /*<! report these most abundant only */
	unsigned int contamination_threshold; /*<! Maximum abundance for contamination */
	int associate_zc; /*<! contamination_threshold = low_bound - 1 */
  unsigned int wsize; /*<! window length */

	/* alignment */
	int score[NUM_NUCLEOTIDES][NUM_NUCLEOTIDES];  /* score matrix for nw alignment */
	int gap_p;  /* penalty for gaps */
	int band;  /* bandwidth for banded nw alignment */
	int ends_free; /* counting the indel at the begining ? Yes[0],No[1] */

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
