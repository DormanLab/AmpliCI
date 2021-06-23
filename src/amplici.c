/**
 * @file amplici.c
 * @author Xiyu Peng
 * 
 * AmpliCI:  Method to identify true sequences (haplotypes) in a sample of
 *           Illumina amplicon data (from FASTQ file).
 * 
 * [TODO] 
 * 1. check reads with unequal lengths and report error
 *    Note the data structure actually supports reads with unequal length
 *    but ampliCI model does not.
 * 2. format your code nicely!  Inspiration?:
 *    https://www.kernel.org/doc/html/v4.10/process/coding-style.html
 * 3. Assigning reads to haplotypes (given or estimated) and output the abundance table 
 * 4. Correct the comments and output in the code. We actually filter on maximum 
 * 		conditional log likelihood (maximum posterior assignment probability), 
 * 		not the reads log likelihood.  
 * 5. use number of observed quality scores, not maxQ- minQ
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <math.h>

#define MATHLIB_STANDALONE 1
// #define DEBUG 0
//#define ADJUST_PVALUE
//#define STORE_FP
#include <Rmath.h>

#include "amplici.h"
#include "ampliclust.h"
#include "initialize.h"
#include "statistics.h"
#include "lmath.h"
#include "io.h"
#include "hash.h"
#include "align.h"
#include "error.h"
#include "error_est.h"

/* amplici */
int amplici_malloc(options *opt, data *dat, initializer *ini, unsigned int **array_fp, double **fp_abun, double **fp_trans, unsigned char **** nw_result, size_t **nw_alen,unsigned int K_space, unsigned n_candidate);
int amplici_realloc(options *opt, initializer *ini, model *mod, unsigned int **array_fp, double **fp_abun, double **fp_trans, unsigned int preK, unsigned int K, unsigned int pre_nfp, unsigned int nfp, size_t sample_size, unsigned int hash_length,unsigned int max_read_length);
void amplici_free(unsigned int *array_fp, double *fp_abun, double *fp_trans, unsigned char *** nw_result, size_t *nw_alen, unsigned int nw_size);

/* estimated scaled true abundance */
int expected_TrueAbundance(options *opt, data *dat, double *H_abun, double *e_trans, double *self_trans,double *abun_true, size_t *idx, unsigned int count_i,unsigned int select, unsigned int i, int conve, double low_bound);
double iterate_expected_true_abundance(unsigned int sample_size, size_t *idx_array, double *e_trans, double *self_trans, double *H_abun, unsigned int select, unsigned int obs_abun, double true_abun);

/* check false positives */
int evaluate_haplotype(options *opt, data *dat, model *mod, initializer *ini, unsigned int K, double low_bound, double * error_profile, unsigned int ord, unsigned int n_candidates, int *fp,int final);
int check_fp_with_indels(options *opt, data *dat, model *mod, initializer *ini,unsigned int select, double low_bound, double * error_profile, int *fp);
int Est_pi(initializer *ini, double *pi, size_t sample_size, unsigned int K, int reassign);
int abun_pvalue(options *opt, initializer *ini, size_t *idx_array, double *e_trans, unsigned int count, unsigned int select, unsigned int threshold, double *p, size_t sample_size, int partial);


/* transition prob with or without alignment */
int Expected_SelfTrans(options *opt, data *dat, double *self_trans,double *error_profile, int err_encoding, double adj);
int ExpTrans_nogap(data *dat, options *opt, initializer *ini,unsigned int H_id, unsigned int select, double *error_profile, int err_encoding);
int ExpTrans_nwalign(data *dat, options *opt, initializer *ini, unsigned char ***nw_result, size_t *nw_alen,unsigned int select, double *error_profile, int err_encoding,unsigned int H_id, double adj);

/**
 * Cluster amplicon sequences.
 * 
 * @param ini	pointer to initializer object
 * @param dat	pointer data object
 * @param opt	pointer to options object
 * @param mod   pointer to model object
 * @param ri    pointer to run_info object
 *
 * return	err status
 **/
int ampliCI(options * opt, data * dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;
	int output_hap = 1;	/* output haplotype FASTA file */
	int use_size = 0;	/* output cluster sizes for abundance */

	/* K to be determined in the function below, opt->K will be changed */
	if ((err = haplotype_selection(opt, dat, mod, ini, opt->K_max))) 
		return err;

	if ((err = realloc_run_info(ri, dat->sample_size, opt->K + 1, 0, 1)))
		return err;
	
	/* assign cluster based on ll */
	assign_clusters(mod->eik, opt->K, dat->sample_size, ri->optimal_cluster_size,
			ri->optimal_cluster_id, 1);

	/* remove reads with too-small log likelihood (Note: ll are alignment-free) */
	/* [KSD] For some reason you do not actually use read log likelihood.
	 * Instead, you use the maximum conditional log likelihood.
	 */
	/* Since it cannot use posterior probabilities. they sum to 1 across all clusters */
	/* [KSD] I was just proposing that this function do as advertised. */
	likelihood_filter(opt->K, opt->ll_cutoff, NULL, mod->pi, ini->e_trans,
							dat->sample_size, ri);

	char *outfile_hap = NULL;
	char *outfile = NULL;
	
	if (opt->outfile_base || opt->outfile_info) {
		FILE *fp = NULL;

		if (!opt->outfile_info) {
			outfile = malloc((strlen(opt->outfile_base) + 5)
							* sizeof (char));
			if (!outfile)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"output file");
   			strcpy(outfile, opt->outfile_base);
			strcat(outfile, ".out");
			opt->outfile_info = outfile;
		}

		fp = fopen(opt->outfile_info, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->outfile_info);

		fprintf(fp, "K: %i\n", opt->K);

		fprintf(fp, "assignments: ");
		fprint_assignment(fp, ri->optimal_cluster_id, dat->sample_size,
								opt->K, 2, 1);
		fprintf(fp, "cluster sizes: ");
		fprint_uints(fp, ri->optimal_cluster_size, opt->K, 3, 1);
		
		fprintf(fp,"pi: ");
		for(unsigned int k = 0; k < opt->K; ++k)
			mod->pi[k] = exp(mod->pi[k]);
		fprint_doubles(fp, mod->pi, opt->K ,6,1);

		fprintf(fp,"reads ll: ");
		fprint_doubles(fp, ri->optimal_cluster_ll, dat->sample_size,
									3, 1);

		//[TODO] output ini->seeds directly
		fprint_fasta(fp, ini->seeds[0], opt->K, 
					 dat->max_read_length, ini->seed_lengths, "H");
		
		fprintf(fp,"ee: ");   // mean expected number of errors
		fprint_doubles(fp, ini->H_ee, opt->K ,3,1);
		
		fprintf(fp, "uniq seq id: ");
		fprint_uints(fp, ini->H, opt->K, 3, 1);

		fprintf(fp, "scaled true abun: ");
		fprint_doubles(fp, ini->H_abun, opt->K, 3, 1);

		fprintf(fp, "obser abun: ");
		for (unsigned k = 0; k < opt->K; k++)
			fprintf(fp, " %*u", 3, ini->uniq_seq_count[ini->H[k]]);
		fprintf(fp, "\n");

		#ifdef ABUN_INTERVAL
		fprintf(fp, "p value: ");
		for (unsigned k = 0; k < opt->K; k++){
			double pvalue = exp(ini->H_pvalue[k]);
			if (pvalue < 1e-3)
				fprintf(fp, " %8.2e", pvalue);
			else
				fprintf(fp, " %.3f", pvalue);
		}
		fprintf(fp, "\n");
		#endif

		if (opt->JC69_model) {
			fprintf(fp, "Estimated common ancestor: \n");
			fprint_fasta(fp, mod->est_ancestor, 1,
					 dat->max_read_length, &dat->max_read_length, "Ancestor");
			fprintf(fp, "Evolution_rate: ");
			fprint_doubles(fp, mod->distance, opt->K, 3, 1);
			fprintf(fp, "log likelihood from JC69 model:%f\n",
								mod->JC_ll);
		}
		fprintf(fp, "log likelihood: %f\n", mod->best_ll);
		fprintf(fp, "Diagnostic Probability threshold: %8.2e\n", opt->p_threshold);
		fprintf(fp, "aic: %f\n", mod->aic);
		fprintf(fp, "bic: %f\n", mod->bic);

		fclose(fp);

		mmessage(INFO_MSG, NO_ERROR, "Output the final result file: "
						"%s \n", opt->outfile_info);
	}
	if(outfile) free(outfile);

	/* format output fasta file for UCHIME */
	if (output_hap) {
		FILE *fp2 = NULL;

		if (!opt->outfile_fasta) {
			if (!opt->outfile_base)
				return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"invalid output filenames");

			outfile_hap = malloc((strlen(opt->outfile_base)
							+ 5) * sizeof(char));
			if (!outfile_hap)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"output file");
 
			strcpy(outfile_hap, opt->outfile_base);
			strcat(outfile_hap, ".fa");
			opt->outfile_fasta = outfile_hap;
		}

		fp2 = fopen(opt->outfile_fasta, "w");
		if (!fp2)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->outfile_fasta);
		
		/* [TODO] rewrote the two function below to allow variable length */
		if (use_size)
			fprint_haplotypes_size(fp2, ini->seeds[0], opt->K,
				dat->max_read_length, opt->p_threshold, "H",
					ini->H_pvalue, ri->optimal_cluster_size,
								ini->H_ee);
		else
			fprint_haplotypes_abun(fp2,ini->seeds[0], opt->K,
				dat->max_read_length, opt->p_threshold, "H", 
					ini->H_pvalue, ini->H_abun, ini->H_ee);

		fclose(fp2);

		mmessage(INFO_MSG, NO_ERROR, "Output the final haplotype fasta "
					"file: %s \n", opt->outfile_fasta);

	}
	
	if(outfile_hap) free(outfile_hap);

	return err;
}/* ampliCI */

/**
 * Select no more than K_max real haplotype sequences 
 * 
 * @param ini	initializer object
 * @param dat	data object
 * @param opt	options object
 * @param mod	model object
 * @param K_max	maximum number of haplotypes selected
 *
 * return	err status
 **/
int haplotype_selection(options * opt, data * dat, model *mod, initializer *ini,
							unsigned int K_max)
{
	int err = NO_ERROR;
	int fxn_debug = opt->info;//DEBUG_I;

	/* ------------------------------------------------------------------ */
	/* Variable Declaration */

	double low_bound; 
	int false_positive;

	double *error_profile = NULL;
	if (opt->use_error_profile && mod->error_profile) {
		error_profile = mod->error_profile;
		debug_msg(DEBUG_II, fxn_debug, "Use error profile. \n");
	}

	if (opt->low_bound > 1) {
		low_bound = opt->low_bound;
	} else {
		low_bound = 2.0;   // avoid singletons 
		mmessage(WARNING_MSG, INVALID_USER_INPUT, "User low bound on "
			"abundance set <= 1: resetting to %f.\n", low_bound);
		if (opt->contamination_threshold != 1) {
			opt->contamination_threshold = 1; 
			mmessage(WARNING_MSG, INVALID_USER_INPUT,
					"resetting contamination threshold %u.",
						opt->contamination_threshold);
		}
	}

	/* malloc space only when we check false positive */
	unsigned int *array_fp = NULL; 
	double *fp_abun = NULL;  // just for test 
	double *fp_trans = NULL;  // store the trans prob 

	unsigned int select = 0;  // num of hap selected
	unsigned int num_fp = 0;  // num of false positive removed 
	unsigned int ord;  // index of current select haplotypes 

	unsigned int K_space = opt->K;       // current space for K clusters 
	unsigned int ini_K = opt->K;        // record the initial K space 
	unsigned int n_candidate = dat->hash_length;

	/* store result from nw alignment between each unique sequences and newly identified sequence */
	unsigned char *** nw_result = NULL;
	size_t *nw_alen = NULL;

	/* Determine number of candidates here, just based on the observed abundance 
	* we will only consider seqs with abundance < low bound */
	
	//double adj_low_bound = low_bound - 0.001;   // try to solve a BUG here
	for (unsigned int i =0; i< dat->hash_length; i++){
		if (low_bound - ini->uniq_seq_count[i] >
			DBL_EPSILON * fmax(ini->uniq_seq_count[i], low_bound)) {    // AVOID A NUMERICAL BUG
			n_candidate = i;
			break;
		}

	}

	if(!n_candidate)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Abundance of all unique sequences are under %f\n",opt->low_bound);

	/* determine the threshold of diagnostic test */
	if (opt->per_candidate)
		opt->p_threshold = opt->alpha / n_candidate;
	else
		opt->p_threshold = opt->alpha;
	
	/* Space mallocation */
	if (((err= amplici_malloc(opt, dat, ini, &array_fp, &fp_abun, &fp_trans,
				&nw_result, &nw_alen, K_space, n_candidate))))
		return err;
	
	/* --------------------------------------------------------------------------- */
	/* Initialization: choose the most abundant unique sequences */

	select = 0;
	num_fp = 0;
	ord = 0;

	/* initialize the abun_true table with observed count */
	for (unsigned int i = 0; i < n_candidate; i++)
		ini->abun_true[i] = ini->uniq_seq_count[i];

	/* calculate the self transition probability */
	Expected_SelfTrans(opt, dat, ini->self_trans, error_profile, mod->err_encoding, mod->adj_trunpois);

	/* select the first haplotype with the highest abundance */
	update_seeds(dat, ini, select, ord);
	
	debug_msg(DEBUG_I, fxn_debug, "Selecting %d with estimated true"
				" abundance %.3f\n", ord, ini->H_abun[select]);

	/* [TODO] parallelize */
	if (opt->nw_align == ALIGNMENT_UNIQ_SEQ)
		ExpTrans_nwalign(dat, opt, ini, nw_result, nw_alen, select,
					error_profile, mod->err_encoding, ord, mod->adj_trunpois);
	else
		ExpTrans_nogap(dat, opt, ini, ord, select, error_profile,
							mod->err_encoding);

	/* use the code only when the haplotype chosen without doubt */
	select++;
	ini->abun_true[ord] = 1.0; // to avoid being selected against

	/* mean expected number of errors for the first haplotype */
	if ((err = mean_exp_errors(dat, ini->uniq_seq_idx[ord],
				ini->uniq_seq_count[ord], &ini->H_ee[0])))
		return err;
	ini->H_pvalue[0] = 0.;  // arbitrarly set to 0

	/* Initialize aic and bic */
	mod->aic = INFINITY;
	mod->bic = INFINITY;
	
	/* seems to work now */
	evaluate_haplotype(opt, dat, mod, ini, 1, low_bound, error_profile,
					0, n_candidate, &false_positive, 1);

	/* ------------------------------------------------------------- */
	/* choose other haplotypes */
	
	while (K_max > 1) {

		false_positive = 0;

		/* If we need more space, we need to realloc the space */
		/* For increasing K */
		if (select == K_space) {
			debug_msg(DEBUG_III, fxn_debug, "begin reallocation");
			K_space = K_space + ini_K;
			if ((err = amplici_realloc(opt, ini, mod, &array_fp,
				&fp_abun, &fp_trans, select, K_space, num_fp,
				num_fp, dat->sample_size, dat->hash_length,
				dat->max_read_length)))
				return err;
			debug_msg(DEBUG_III, fxn_debug,
						"Finish reallocation\n");
		}

		/* For increasing num_fp */
		#ifdef STORE_FP
		if (opt->check_false_positive && num_fp == fp_space) {
			fp_space = fp_space * 2;
			if ((err = amplici_realloc(opt, ini, mod, &array_fp,
				&fp_abun, &fp_trans, K_space, K_space, num_fp,
				fp_space, dat->sample_size, dat->hash_length,
							dat->max_read_length)))
				return err;
		}
		#endif

		/* update relative true abundance and 
		 * choose a haplotype with the highest relative true abundance
		 */
		double max = 0.;
		/* [TODO] parallelize */
		//#pragma omp parallel for
		for (unsigned int i = 0; i < n_candidate; i++) {

			if (ini->abun_true[i] >= low_bound) {

				if ((err =  expected_TrueAbundance(opt, dat,
					ini->H_abun, ini->e_trans,
					ini->self_trans, ini->abun_true,
					ini->uniq_seq_idx,
					ini->uniq_seq_count[i], select, i,
					opt->convergence_amplici, low_bound)))
					mmessage(WARNING_MSG, INTERNAL_ERROR,
						"warning about no convergence "
						"generated when updating true "
						"abundance of %i in %ith step \n",
								i, select + 1);

				debug_msg(DEBUG_III, fxn_debug, "Estimated "
					"abundance of %d: %.3f\n", i,
					ini->abun_true[i]);
			}
		}
		
		for (unsigned int i =0; i < n_candidate; i++)
			if (ini->abun_true[i] > max) {
				ord = i;
				max = ini->abun_true[i];
			}

		if (ini->abun_true[ord] < low_bound)
			break;

		/* select the haplotype temporarily */
		update_seeds(dat, ini, select, ord);
	
		debug_msg(DEBUG_I, fxn_debug, "Selecting %d with estimated true"
				" abundance %.3f\n", ord, ini->H_abun[select]);
			
		/* Transition prob without alignment free strategy */
		if (opt->nw_align == ALIGNMENT_UNIQ_SEQ)
			ExpTrans_nwalign(dat, opt, ini, nw_result, nw_alen,
				select, error_profile, mod->err_encoding, ord, mod->adj_trunpois);
		else if (opt->nw_align == NO_ALIGNMENT)
			ExpTrans_nogap(dat, opt, ini, ord, select,
					error_profile, mod->err_encoding);
		/* computed in evaluate_haplotype() for opt->nw_align ==
		 * ALIGNMENT_HAPLOTYPES
		 */
		

		/* evaluate the newly selected haplotype 
		 * check if it is a false positive
		 */
		if (opt->check_false_positive) {
			if ((err = evaluate_haplotype(opt, dat, mod, ini, 
				select + 1, low_bound, error_profile, ord,
					n_candidate, &false_positive, 0)))
				return err;

			#if DEBUG
			debug_msg(DEBUG_I, fxn_debug, "check the %d th unique sequence\n", ord);
			for (size_t j = 0; j < dat->max_read_length; ++j)
				fprintf(stderr, "%c", xy_to_char[ini->seeds[select][j]]);

			fprintf(stderr, "\n");
			#endif

			if (false_positive) {
				#ifdef STORE_FP
				array_fp[num_fp] = ord;
				fp_abun[num_fp] = ini->abun_true[ord];
				double * des_sta = &fp_trans[num_fp * dat->sample_size];
				double * ori_sta = &ini->e_trans[select * dat->sample_size];
				memcpy(des_sta, ori_sta, dat->sample_size * sizeof *des_sta);
				#endif

				num_fp ++;
				ini->abun_true[ord] = 1.0;  // To avoid being selected again 

				debug_msg(DEBUG_I, fxn_debug, "remove %d since "
						"it is a false positive\n", ord);
				continue;
			}
		}

		select ++;
		ini->abun_true[ord] = 1.0;

		if (select == K_max)
			break;

	}

	opt->K = select;

	// update E matrix, mod->ll, mod->param, mod->bic (modified) , mod->aic (modified)
	debug_msg(DEBUG_I, fxn_debug, "Final Statistics: \n");
	evaluate_haplotype(opt, dat, mod, ini, opt->K, low_bound, error_profile,
					0, n_candidate, &false_positive, 1);

	for (unsigned int k = 0; k < opt->K; k++) {

		debug_msg(SILENT, fxn_debug, "True abundance of the %ith "
			"haplotype (%ith unique seq id): %.3f\n", k, ini->H[k],
			ini->H_abun[k]);
		debug_msg(SILENT, fxn_debug, "observed abundance of the %ith "
			"haplotype (%i th unique seq id): %i\n", k, ini->H[k],
			ini->uniq_seq_count[ini->H[k]]);
		debug_msg(DEBUG_I, fxn_debug, "The mean exp errors of the %ith "
			"haplotype (%i th unique seq id): %.3f\n", k, ini->H[k],
			ini->H_ee[k]);
	}

	mod->best_ll = mod->ll;

	amplici_free(array_fp, fp_abun, fp_trans, nw_result, nw_alen,
							dat->hash_length);
	opt->K_space = K_space;

	return err;
}/* haplotype_selection */

/**
 * Free space if we store arrays for false positives or use ALIGNMENT_UNIQ_SEQ
 * 
 * @param array_fp  sequences of false positives
 * @param fp_abun 	estiamted abundance of false positives
 * @param fp_trans  transition prob of false positive
 * @param nw_result alignment result
 * @param nw_alen 	alignment length
 * @param nw_size   num. of alignment result 
 *
 **/
void amplici_free(unsigned int *array_fp, double *fp_abun, double *fp_trans,
					unsigned char *** nw_result, size_t *nw_alen, unsigned int nw_size){
	if(array_fp)free(array_fp);
	if(fp_abun)free(fp_abun);
	if(fp_trans)free(fp_trans);
	if (nw_result){
		free_nw_result(nw_result, nw_size);
		free(nw_result);
	}
	if(nw_alen)free(nw_alen);
	
}/* amplici_free */

/**
 * Free space used for NW alignments.
 * 
 * @param nw_result	the memory with alignments
 * @param space		number of alignments
 *
 **/
void free_nw_result(unsigned char ***nw_result, unsigned int space)
{
	for (unsigned int i = 0; i < space; i++) {
		if (nw_result[i]) {
			for (unsigned int j = 0; j < 2; j++)
				if (nw_result[i][j])
					free(nw_result[i][j]);
			free(nw_result[i]);
		}
		nw_result[i] = NULL;
	}
	//free(nw_result);
}/* free_nw_result */


/**
 * Obtain Needleman-Wunsch alignment between a selected haplotype and all other
 * unique sequences.  Stores alignment for each unique sequence in
 * \par nw_result and alignment length in \par nw_alen.  Also stores number
 * of mismatches in initializer::nw_mismatch and number of indels in
 * initializer::nw_indels.
 *
 * [KSD] This code is largely duplicated with check_fp_with_indels().  Consider
 * [KSD] putting the duplicated stuff in one location.  Also, it is not
 * [KSD] clear why the results are stored in the initializer.  In fact, it is
 * [KSD] not clear why there is an initializer object at all.  Put in the data
 * [KSD] object?
 * 
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param ini		pointer to initializer object
 * @param nw_result	nw alignment result (to be set)
 * @param nw_alen	length of nw alignemnt (to be set)
 * @param size		num of unique sequences
 * @param select	current selected haplotype
 * 
 * @return err		error status
 **/
int nwalign_matrix(options *opt, data *dat, initializer *ini,
	unsigned char ***nw_result, size_t *nw_alen, unsigned int size,
							unsigned int select)
{

	int err = NO_ERROR;
	//unsigned int size;
	size = dat->hash_length;
	int ends_free = 1;  // should be ends_free alignment. No panalty for gaps in the end 
	
	unsigned char *ref_seq = ini->seeds[select];
	unsigned int ref_len = ini->seed_lengths[select];

	/* align each unique sequence to current haplotype */
	for (unsigned int i = 0; i < size; ++i) {
		size_t alen = dat->max_read_length;
		unsigned char **aln = NULL;
		unsigned int rlen = dat->lengths[i];

		/* only when sequences are different, we need alignment */
		aln = nwalign(ref_seq, dat->dmat[ini->uniq_seq_idx[i]],
			(size_t) ref_len,
			(size_t) dat->lengths[ini->uniq_seq_idx[i]],
			opt->score, opt->gap_p, opt->band, ends_free, NULL,
								&err, &alen);
		
		#if DEBUG
		fprintf(stderr, "alignment");
		for (size_t j = 0; j < alen; ++j) {
				fprintf(stderr, "%c", aln[0][j] == '-'
				? '-' : xy_to_char[(int) aln[0][j]]);
			}
			fprintf(stderr, "\n");
			for (size_t j = 0; j < alen; ++j) {
				fprintf(stderr, "%c", aln[1][j] == '-'
				? '-' : xy_to_char[(int) aln[1][j]]);
			}
		#endif
		
		nw_result[i] = aln;
		nw_alen[i] = alen;

		/* calculate number of indels and mismatch based on alignment */
		ana_alignment(aln, alen, rlen, &ini->nw_indels[select*size+i], 
					&ini->nw_mismatch[select*size+i], opt->ends_free,opt->info);
	}

	return err;
}/* nwalign_matrix */

/**
 * malloc additional space specific for amplici 
 * 
 * @param opt   	pointer to options object
 * @param dat 		pointer to data object
 * @param ini   		pointer to initializer object
 * @param array_fp  	sequences of false positives
 * @param fp_abun 		estiamted abundance of false positives
 * @param fp_trans  	transition prob of false positive
 * @param nw_result 	nw alignment result 
 * @param nw_alen   	length of nw alignemnt
 * @param K_space   	space for K clusters
 * @param n_candidate   total number of candidates
 * 
 * @return err      error status
 **/
int amplici_malloc(options *opt, data *dat, initializer *ini,
	unsigned int **array_fp, double **fp_abun, double **fp_trans,
	unsigned char ****nw_result, size_t **nw_alen, unsigned int K_space,
							unsigned n_candidate)
{	
	UNUSED(array_fp);
	UNUSED(fp_abun);
	UNUSED(fp_trans);

	/* true abundance table */
	if (!ini->abun_true)
	    	ini->abun_true = malloc(n_candidate * sizeof * ini->abun_true);
	if (!ini->abun_true)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.abun_true");

	/* haplotypes idx in unique sequence table 
	 * need reallocation when K increases
	 */
   	if (!ini->H)
   		ini->H = malloc(K_space * sizeof * ini->H);
	if (!ini->H)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.H");

	/* Haplotypes abundance */
	if (!ini->H_abun)
		ini->H_abun = malloc(K_space * sizeof * ini->H_abun);
	if (!ini->H_abun)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.H_abun");

	/* haplotypes expected number of errors */
	if (!ini->H_ee)
		ini->H_ee = malloc(K_space * sizeof * ini->H_ee);
	if (!ini->H_ee)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.H_ee");
	
	/* p-value of haplotypes */
	if (!ini->H_pvalue)
		ini->H_pvalue = malloc(K_space * sizeof *ini->H_pvalue);
	if (!ini->H_pvalue)
		return mmessage(ERROR_MSG,MEMORY_ALLOCATION,"amplici.H_pvalue");

	/* false positive table may need reallocation */
#ifdef STORE_FP
	if (opt->check_false_positive) {
		//if (!*array_fp)
    		*array_fp = malloc(K_space * sizeof **array_fp);
		//if (!*fp_abun)
		*fp_abun = malloc(K_space * sizeof **fp_abun);
		//if (!*fp_trans)
		*fp_trans = malloc(K_space * dat->sample_size * sizeof **fp_trans);

		if (!*fp_abun || !*array_fp || !*fp_trans )
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.fp");
	}
#endif
	
	/* malloc transition matrix between reads and haplotypes (log value) 
	 * need reallocation when K increases
	 */
	if (!ini->e_trans)
		ini->e_trans = malloc(dat->sample_size * K_space * sizeof * ini->e_trans);
	if (!ini->e_trans)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.e_trans");

	/* self transition matrix */
	if (!ini->self_trans)
   		ini->self_trans = malloc(dat->sample_size * sizeof * ini->self_trans);
	if (!ini->self_trans)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.self_trans");

	/* simplified results of nw alignment, (record it initializer struct )
	need reallocation when K increases.
	could also be used as a bound to update true abundance. */
	if (opt->nw_align == ALIGNMENT_UNIQ_SEQ) {
		if (!ini->nw_mismatch)
			ini->nw_mismatch = calloc(dat->hash_length * K_space, sizeof * ini->nw_mismatch);
		if (!ini->nw_mismatch)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.nw_mismatch");
		if (!ini->nw_indels)
			ini->nw_indels = calloc(dat->hash_length * K_space, sizeof * ini->nw_indels);
		if (!ini->nw_indels)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.nw_indels");

	/* results of nw alignment */
	/* should think about a more effcient way to do this. nw alignment is really time consuming */
	/* Just store the result for each haplotypes (memory efficient) */
		//if (!*nw_result)
			*nw_result = malloc(dat->hash_length * sizeof(**nw_result));
		//if (!*nw_alen)
			*nw_alen = malloc(dat->hash_length * sizeof ** nw_alen);

		if (!*nw_result || !*nw_alen)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "amplici.nw");
	}
	return NO_ERROR;
}/* amplici_malloc */

/**
 * realloc space when num of clusters or num of false positives increases
 * 
 * @param opt			pointer to options object
 * @param mod			pointer to model object
 * @param ini 			pointer to initializer object
 * @param array_fp		sequences of false positives
 * @param fp_abun		estimated abundance of false positives
 * @param fp_trans		transition prob of false positive
 * @param preK			current space for preK clusters
 * @param K			realloc space for K clusters
 * @param pre_nfp		current space for pre_nfp false positives
 * @param nfp			realloc space for nfp false positives 
 * @param sample_size		number of total reads
 * @param hash_length		total number of candidates
 * @param max_read_length	length of reads
 * 
 * @return			error status
 * 
 * need to realloc space for pi, haplotype, distance, eik in model object
 * seeds, seed_idx, seed_lengths in initializer object
 * */
int amplici_realloc(options *opt, initializer *ini, model *mod,
	unsigned int **array_fp, double **fp_abun, double **fp_trans,
	unsigned int preK, unsigned int K, unsigned int pre_nfp,
	unsigned int nfp, size_t sample_size, unsigned int hash_length,
	unsigned int max_read_length)
{	
	int err = NO_ERROR;


#ifdef STORE_FP
	unsigned int fp_change = nfp - pre_nfp;
	if (fp_change) {
		unsigned int *array_fp_copy = realloc(*array_fp, nfp * sizeof **array_fp);
		double * fp_abun_copy = realloc(*fp_abun, nfp * sizeof **fp_abun);
		double *fp_trans_copy = realloc(*fp_trans,nfp * sample_size * sizeof **fp_trans);
		if (!array_fp_copy || !fp_abun_copy || !fp_trans_copy) {
			if (array_fp_copy) free(array_fp_copy);
			if (fp_abun_copy) free(fp_abun_copy);
			if (fp_trans_copy) free(fp_trans_copy);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"amplici.realloc.fp");
		}
		*array_fp = array_fp_copy;
		*fp_abun = fp_abun_copy;
		*fp_trans = fp_trans_copy;
	}
#else
	UNUSED(array_fp);
	UNUSED(fp_abun);
	UNUSED(fp_trans);
	UNUSED(nfp);
	UNUSED(pre_nfp);
#endif

	unsigned int K_change = K - preK;

	if (K_change > 0) {

		/* H and H_abun */
		unsigned int *H = realloc(ini->H, K * sizeof * ini->H);
		double *H_abun = realloc(ini->H_abun, K * sizeof * ini->H_abun);
		double *H_ee = realloc(ini->H_ee, K * sizeof * ini->H_ee);
		double *H_pvalue = realloc(ini->H_pvalue,K * sizeof * ini->H_pvalue);
		
		if (!H || !H_abun || !H_ee || !H_pvalue) {
			if (H) free(H);
			if (H_abun) free(H_abun);
			if (H_ee) free(H_ee);
			if (H_pvalue) free(H_pvalue);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"amplici.realloc.H");
		}
		ini->H = H;
		ini->H_abun = H_abun;
		ini->H_ee = H_ee;
		ini->H_pvalue = H_pvalue;

		/* e_trans */
		double *e_trans = realloc(ini->e_trans, sample_size * K * sizeof *ini->e_trans);
		if (!e_trans)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"amplici.realloc.e_trans");
		ini->e_trans = e_trans;

		/* nw mismatch and mw_indels */
		if (opt->nw_align == ALIGNMENT_UNIQ_SEQ) {
			unsigned int *nw_mismatch = realloc(ini->nw_mismatch,
				hash_length * K * sizeof *ini->nw_mismatch);
			unsigned int *nw_indels = realloc(ini->nw_indels,
				hash_length * K * sizeof * ini->nw_indels);
			if (!nw_mismatch || !nw_indels) {
				if (nw_mismatch) free(nw_mismatch);
				if (nw_indels) free(nw_indels);
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"amplici.realloc.nw");
			}
			ini->nw_mismatch = nw_mismatch;
			ini->nw_indels = nw_indels;
		}

		/* haplotype, pi, eik */
		double *pi = realloc(mod->pi, K * sizeof *mod->pi);
		//unsigned char *haplotypes = realloc(mod->haplotypes,
		//	max_read_length * K * sizeof *mod->haplotypes);  
		double *eik = realloc(mod->eik,
					sample_size * K * sizeof *mod->eik);
		double *distance = realloc(mod->distance,
						K * sizeof *mod->distance);
		double *JC_ll_K = realloc(mod->JC_ll_K,
						K * sizeof *mod->JC_ll_K);
 		unsigned int *cluster_size = realloc(ini->cluster_size,
						K * sizeof *ini->cluster_size);

		if (!pi || !eik || !distance || !cluster_size
								|| !JC_ll_K) {
			if (pi) free(pi);
			//if (haplotypes) free(haplotypes);
			if (eik) free(eik);
			if (distance) free(distance);
			if (cluster_size) free(cluster_size);
			if (JC_ll_K) free(JC_ll_K);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"amplici.realloc.mod");
		}
		mod->pi = pi;
		//mod->haplotypes = haplotypes;
		mod->eik = eik;
		mod->distance = distance;
		mod->JC_ll_K  = JC_ll_K;
		ini->cluster_size = cluster_size;

		if((err = realloc_seeds(ini, max_read_length, preK, K)))
			return err;

		/* seeds, seeds_length, seed_idx */
		//size_t *seed_idx = realloc(ini->seed_idx, K * sizeof *ini->seed_idx);   
		/* 
		unsigned int *seed_lengths = realloc(ini->seed_lengths, K * sizeof *ini->seed_lengths);
		data_t **seeds = realloc(ini->seeds, K * sizeof *ini->seeds); 

		if ( !seeds || !seed_lengths) {
			//if (seed_idx) free(seed_idx);
			if (seeds) free(seeds);
			if (seed_lengths) free(seed_lengths);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"amplici.realloc.seed");
		}
		// ini->seed_idx = seed_idx;
		ini->seed_lengths = seed_lengths;   // Uninitialized 
		ini->seeds = seeds;

		data_t *dptr = realloc(ini->seeds[0], max_read_length * K * sizeof **ini->seeds);
		if (!dptr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"reallloc.initializer.seeds");
		size_t s = 0;
		if (ini->seeds[0] == dptr)  s = preK;

		for (size_t k = s; k < K; k++) {
			ini->seeds[k] = dptr;
			dptr += max_read_length;
		} */
		/*
		for (unsigned int k = preK; k < K; k++) {
			seeds[k] = NULL;
			seeds[k] = calloc(max_read_length, sizeof **ini->seeds);    
			if (!seeds[k]) {
				//free(seed_idx);
				free(seeds);
				free(seed_lengths);
				for (unsigned int i = preK; i < k; ++i)
					free(seeds[i]);
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"amplici.realloc.seeds");
			}
		}
		ini->seeds = seeds;
		*/
	}

	return err;

}/* amplici_realloc */

/**
 * Update estimated scaled true abundance for each uniq seq
 *
 * @param opt		pointer to opt object
 * @param dat		pointer to data object
 * @param H_abun	abundance of haplotypes
 * @param self_trans	log probability of misread to each read from itself
 * @param abun_true	estimated true abundances of all uniq seq
 * @param e_trans	log probability of misread to each read from existing
 *			haplotypes
 * @param idx		idx of the unique seq in dmat 
 * @param i		idx of current unique sequence 
 * @param select	total number of selected haplotypes
 * @param count_i	observed abundance of unique sequence
 * @param conve		convergence or not when updating abundance
 * @param low_bound	lower bound of haplotype abundance
 *
 * @return		error status
 */
int expected_TrueAbundance(options *opt, data *dat, double *H_abun,
	double *e_trans, double *self_trans, double *abun_true, size_t *idx,
	unsigned int count_i, unsigned int select, unsigned int i, int conve,
							double low_bound)
{

	int err = NO_ERROR;

	double delta;
	unsigned int iter = 0;
	
	/* find the read idx in the sample set of this unique sequence */
	//hash *unit = NULL;
	unsigned char *seq = dat->dmat[idx[i]];
	unsigned int length = dat->lengths[idx[i]];
	//HASH_FIND(hh, dat->seq_count, seq, length * sizeof *seq, unit);
	size_t *idx_array = NULL;
	if (((err = find_index(dat->seq_count, seq, length, &idx_array))))
		return err;

	/* update true abundance of unique sequence i */
	double true_abun_i = iterate_expected_true_abundance(dat->sample_size,
					idx_array, e_trans, self_trans, H_abun,
						select, count_i, abun_true[i]);
	iter++;

	if (conve) {
		do {
			delta = (abun_true[i] - true_abun_i) / abun_true[i];
			//delta = abun_true[i]-true_abun_i;
			/* updated abundance should be smaller than current abundance */
			
			// previous if((delta < 0) || (true_abun_i < 0))
			// have fixed the bug in numerical calculations. 
			abun_true[i] = true_abun_i;

			if (fabs(delta) < opt->epsilon
						|| true_abun_i < low_bound)
				return NO_ERROR;

			if (delta < 0)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
						"true abundance increase.\n");

			true_abun_i = iterate_expected_true_abundance(
					dat->sample_size, idx_array,
					e_trans, self_trans, H_abun, select,
							count_i, abun_true[i]);
			if (true_abun_i < 0)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
						"True abundance under 0.\n");

			iter ++;

			if (iter > opt->n_iter_amplici)
				return mmessage(WARNING_MSG, INTERNAL_ERROR,
					"exceed max interations,%i\n", iter);

		} while(1);
	} else {
		abun_true[i] = true_abun_i;
	}

	return NO_ERROR;
}/* expected_TrueAbundance */

/**
 * One iterate of the fixed point iteration to compute expected abundance of
 * a unique sequence.
 *
 * @param sample_size	number of total reads
 * @param idx_array  	idx of all reads of the uniq seq in dmat 
 * @param e_trans	log probability of misread to each read from existing
 *			haplotypes
 * @param self_trans	log probability of misread to each read from itself
 * @param H_abun	abundance of haplotypes
 * @param select	total number of selected haplotypes
 * @param obs_abun	observed abundance of unique sequence
 * @param true_abun	estimated expected true abundance of unique sequence, x
 * 
 * @return		new estimated expected true abundance, f(x)
 */
double iterate_expected_true_abundance(unsigned int sample_size, size_t *idx_array,
	double *e_trans, double *self_trans, double *H_abun, unsigned int select,
	unsigned int obs_abun, double true_abun){

	double sum_pro_H;
	double true_abun_new  = obs_abun;

	for (unsigned int r = 0; r < obs_abun; r++) {
		/* prob rmi <- H */
		sum_pro_H = 0.;
		for (unsigned int k =0 ; k < select; k++)
			sum_pro_H += exp(e_trans[k * sample_size
					+ idx_array[r]]) * H_abun[k];
		true_abun_new -= sum_pro_H / (sum_pro_H
			+ exp(self_trans[idx_array[r]]) * true_abun);
	}
	
	return true_abun_new;
} /* iterate_expected_true_abundance */

/**
 * Compute transition probability between all reads and current haplotype
 * based on nw alignment result.
 * 
 * @param dat			pointer to data object
 * @param opt			pointer to options object
 * @param ini 			pointer to initializer object
 * @param nw_result		nw alignment result
 * @param nw_alen		length of nw alignment
 * @param select		current number of haplotypes
 * @param error_profile		input error profile
 * @param error_encoding	nucleotide encoding of error profile
 * @param H_id			idx of the current haplotype candidate
 * @param adj			Pr(#{indels} < read_length)
 * 
 * @return err			error status
 **/
int ExpTrans_nwalign(data *dat, options *opt, initializer *ini,
	unsigned char ***nw_result, size_t *nw_alen, unsigned int select,
	double *error_profile, int err_encoding, unsigned int H_id, double adj)
{

	int err = NO_ERROR;

	/* align every unique sequence to the candidate haplotype */
	if ((err = nwalign_matrix(opt, dat, ini, nw_result, nw_alen,
						dat->hash_length, select)))
		return err;

	unsigned int n = dat->sample_size;
	double *e_trans = ini->e_trans;
	
	if (!e_trans)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"cal_e_trans_nw.e_trans");

	for (unsigned int r = 0; r < n; r++) {

		unsigned int id = ini->reads_uniq_id[r];   // id in unique sequenes table 
		unsigned int idx = n * select + r;

		/* transition prob between same sequences with different 
		 * quality score; no need for alignment in this case
		 */
		if (id == H_id) {
			e_trans[idx] = ini->self_trans[r];  // equal to self transition probability
		} else {

			unsigned char **align = nw_result[id];
			size_t alen = nw_alen[id];

			e_trans[idx] = trans_nw(opt, align, alen,
				ini->nw_mismatch[dat->hash_length*select+id],
				ini->nw_indels[dat->hash_length*select+id],
				error_profile, err_encoding, dat->qmat[r],
					dat->n_quality, adj, dat->lengths[r],
							dat->error_prob,opt->ends_free);
		}
	}

	/* free the alignments */
	// realloc_nw_result(nw_result, dat->hash_length);
	free_nw_result(nw_result, dat->hash_length); // [KSD] recommends

	return err;
}/* ExpTrans_nwalign */

/**
 * Self transition probability 
 * 
 * @param opt			pointer to options object
 * @param dat			pointer to data object
 * @param self_trans		self transition prob array
 * @param error_profile		input error profile
 * @param error_encoding	nucleotide encodings of error profile
 * @param adj			Pr(#{indels} <= read_length)
 * @return err			error status
 **/
int Expected_SelfTrans(options *opt, data *dat, double *self_trans,
			double *error_profile, int err_encoding, double adj)
{

	double log_epsilon = opt->epsilon_aln;

	/* [TODO] parallelize */
	//#pragma omp parallel for
	for (unsigned int r = 0; r < dat->sample_size; r++) {
		self_trans[r] = 0;		
			
		if (opt->indel_model == INDEL_PER_READ
					&& opt->nw_align == ALIGNMENT_UNIQ_SEQ)
			self_trans[r] += dpois(0, opt->indel_error
						* dat->lengths[r], 1) - adj;

		for (unsigned int j = 0; j < dat->lengths[r]; j++) {
			
			/* to accelerate? */
//			if( self_trans[r] < log_epsilon ) {
//				self_trans[r] = log_epsilon;
//				break;
//			}
		
			if (error_profile) {
				if (err_encoding == STD_ENCODING)
					self_trans[r] +=
						translate_error_STD_to_XY(error_profile,
						dat->n_quality, dat->dmat[r][j],
						dat->dmat[r][j], dat->qmat[r][j]);
				else if (err_encoding == XY_ENCODING)
					self_trans[r] += 
						error_profile[(NUM_NUCLEOTIDES
						* dat->dmat[r][j] + dat->dmat[r][j])
						* dat->n_quality+dat->qmat[r][j]];
			} else {
				//double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);	
				double ep = dat->error_prob[dat->qmat[r][j]];					
				self_trans[r] += log(1 - ep);
			}
		}
	}
	return NO_ERROR;
}/* Expected_SelfTrans */

/**
 * Transition probability without nw alignment.
 * 
 * @param opt   		pointer to options object
 * @param dat 			pointer to data object
 * @param ini   		pointer to initializer object
 * @param select		total number of selected haplotypes
 * @param H_id 			idx of the current haplotype candidate
 * @param error_profile  	input error profile
 * @param error_encoding 	nucleotide encodings of error profile
 * 
 * @return err      		error status
 **/
int ExpTrans_nogap(data *dat, options *opt, initializer *ini, unsigned int H_id,
		unsigned int select, double *error_profile, int err_encoding)
{

//	double log_epsilon = opt->epsilon_aln;
	unsigned char *seq = ini->seeds[select];

	unsigned int n = dat->sample_size;
	unsigned int idx;
	double l1third =  log(1./3);
				
	for (unsigned int r = 0; r < dat->sample_size; r++) {
		unsigned int id = ini->reads_uniq_id[r]; 
		idx = n * select + r;

		if (id == H_id) {
			ini->e_trans[idx] = ini->self_trans[r]; 
		} else {

			ini->e_trans[idx] = 0.; 
			for (unsigned int j = 0; j < dat->lengths[r]; j++) {
			
				/* to accelerate? */
//				if (ini->e_trans[idx] < log_epsilon) {
//					ini->e_trans[idx] = log_epsilon;
//					break;
//				}

				if (error_profile) {
					if (err_encoding == STD_ENCODING)
						ini->e_trans[idx] +=
							translate_error_STD_to_XY(
							error_profile,
							dat->n_quality, seq[j],
							dat->dmat[r][j],
							dat->qmat[r][j]);
					else if (err_encoding == XY_ENCODING)
						ini->e_trans[idx] += 
							error_profile[(NUM_NUCLEOTIDES
							* seq[j] + dat->dmat[r][j])
							* dat->n_quality
								+ dat->qmat[r][j]];
				} else {
					//double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);
					double ep = dat->error_prob[dat->qmat[r][j]];	
					if (dat->dmat[r][j] == seq[j])			
						ini->e_trans[idx] += log(1 - ep);
					else
						ini->e_trans[idx] += log(ep)
								+ l1third;
				}
			}
		}
	}
	return NO_ERROR;
}/* ExpTrans_nogap */


/**
 * Transition probability (logged) between one read and current haplotype 
 * based on nw alignment result.
 * 
 * @param opt			pointer to options object
 * @param aln			nw alignment result
 * @param alen 			alignment length
 * @param mismatch		num of mismatches ([KSD] Why mismatch?)
 * @param ngap 			num of gaps 
 * @param error_profile		input error profile
 * @param error_encoding	nucleotide encodings of error profile
 * @param rqmat			read quality score sequence
 * @param adj			Pr(#{indels} <= read_length)
 * @param rlen			read length 
 * @param n_quality		the number of possible quality scores
 * @param error_prob	error prob indicated by quality score
 * @param ends_free     count the indel at the beginning of the alignment ?
 * 						Yes[0], No[1]
 * 
 * @return e_trans		log transition prob
 **/
double trans_nw(options *opt, unsigned char **aln, size_t alen,
	unsigned int mismatch, unsigned int ngap, double *error_profile,
	int err_encoding, unsigned char *rqmat, unsigned char n_quality, 
	double adj, unsigned int rlen, double *error_prob, int ends_free)
{
		
	int fxn_debug = ABSOLUTE_SILENCE;
	double log_epsilon = opt->epsilon_aln;
	double e_trans = 0;

	/* prob of an indel event */
	if (opt->indel_model == INDEL_PER_READ)
		e_trans += dpois(ngap, opt->indel_error * rlen, 1) - adj;

	/* may be used later */
	unsigned int pre_nmismatch = mismatch;
	unsigned int pre_ngap = ngap;

	double l1third = log(1./3);

	unsigned int nins = 0; 
	unsigned int ndel = 0;
	unsigned int nindel = 0;
	unsigned int nmatch = 0;
	unsigned int nmismat = 0;
	
	if (aln) {
		for (size_t j = 0; j < alen; j++) {
		
// commented Nov 25
//			if (e_trans < log_epsilon) {   // for acceleration 
//				e_trans = log_epsilon;
//				break;
//			}

			unsigned int j1 = j - nins;   // pos idx of hap
			unsigned int j2 = j - ndel;   // pos idx of read

			/* gaps in the end */
			if (j2 >= rlen)   // No data
				break;

			if (j1 >= rlen) {    //assume there is no error for the missing data in the end
				if (opt->indel_model != INDEL_PER_SITE1
					&& opt->indel_model != INDEL_PER_SITE2)
					continue;

				/* [KSD] The gory details of computing a
				 * [KSD] read error should not be
				 * [KSD] exposed here and other places.
				 * [KSD] Consider an (inline) function.
				 */

				/* for per-site indel models */
				if (error_profile) {
					if (err_encoding == STD_ENCODING)
						e_trans += translate_error_STD_to_XY(
							error_profile,
							n_quality, aln[1][j],
							aln[1][j], rqmat[j2]);
					else if (err_encoding == XY_ENCODING)
						e_trans += error_profile[
							(NUM_NUCLEOTIDES
							* aln[1][j]
							+ aln[1][j])
							* n_quality
							+ rqmat[j2]];
				} else {
					double ep = error_prob[rqmat[j2]];
					e_trans += log(1 - ep);
				}
				continue;
			}

			/* insertion */
			if (aln[0][j] == '-') {
				nins ++;

				if (j == 0)
					nindel += ends_free ? 0: 1;
				else if (aln[0][j-1] != '-')
					nindel++;

				if (opt->indel_model == INDEL_PER_SITE1)
					e_trans += log(opt->insertion_error);	/* [KSD] precompute */
				else if (opt->indel_model == INDEL_PER_SITE2)
					e_trans += log(1./4);	/* [KSD] precompute */
				continue;
			}

			/* deletion */
			if (aln[1][j] == '-') {
				ndel ++;
				if (j == 0)
					nindel += ends_free ? 0: 1;
				else if (aln[1][j-1] != '-')
					nindel++;
				if (opt->indel_model == INDEL_PER_SITE1)
					e_trans += log(opt->deletion_error);
				continue;
			}

			/* match and mismatch */
			if (error_profile) {
				if (err_encoding == STD_ENCODING)
					e_trans += translate_error_STD_to_XY(
						error_profile, n_quality,
						aln[0][j], aln[1][j],
							rqmat[j2]);
				else if (err_encoding == XY_ENCODING)
					e_trans += error_profile[(
						NUM_NUCLEOTIDES
						* aln[0][j] + aln[1][j])
						* n_quality + rqmat[j2]];
			
				if (aln[0][j] == aln[1][j])
					nmatch++;
				else
					nmismat++; 
			} else {
				double ep = error_prob[rqmat[j2]];

				if (aln[0][j] == aln[1][j]) {
					nmatch++;
					e_trans += log(1 - ep);
				} else {
					nmismat++;
					e_trans += log(ep) + l1third;
				}
			}
		}
	/* [KSD] This is when banded alignment fails? YES */
	} else {
		e_trans = log_epsilon;
	}

	/* codes for debugging */
	if (fxn_debug & pre_ngap & pre_nmismatch) {
		if (nindel != pre_ngap || nmismat != pre_nmismatch) {
			mmessage(INFO_MSG, NO_ERROR, "alen: %d\n", alen);
			mmessage(INFO_MSG, NO_ERROR, "pre_ngap: %d; pre_mismatch:%d \n", pre_ngap,pre_nmismatch);
			mmessage(INFO_MSG, NO_ERROR, "ngap: %d; mismatch: %d \n", nindel,nmismat);
			mmessage(ERROR_MSG, INTERNAL_ERROR, "result not match, should have bug here\n");
			for (size_t j = 0; j < alen; ++j)
				fprintf(stderr, "%c", aln[0][j] == '-'
					? '-' : xy_to_char[(int) aln[0][j]]);
			fprintf(stderr, "\n");
			for (size_t j = 0; j < alen; ++j)
				fprintf(stderr, "%c", aln[1][j] == '-'
					? '-' : xy_to_char[(int) aln[1][j]]);
		}
	}

	return e_trans;
}/* trans_nw */

/**
 * Decide whether a candidate haplotype should be included in the haplotype set.
 * - If \par options::nw_align == ALIGNMENT_HAPLOTYPES, then align the candidate
 *   to existing haplotypes and
 *   Check false positives or compute model log likelihoood and bic (final)
 * 
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param K		number of selected haplotypes, including candidate
 * @param low_bound	lower bound of haplotype abundance
 * @param error_profile	input error profile
 * @param ord		idx of the current haplotype candidate
 * @param n_candidates	number of total candidates
 * @param fp		false positive or not (return)
 * @param final		final step ?
 * 
 * @return		error status
 * 
 * * [TODO] it is too complex. Simplify it 
 **/
int evaluate_haplotype(options *opt, data *dat, model *mod, initializer *ini,
	unsigned int K, double low_bound, double *error_profile,
	unsigned int ord, unsigned int n_candidates, int *fp, int final)
{
		
	int fxn_debug = opt->info;
	int err = NO_ERROR;
	double new_bic = INFINITY, pre_bic = mod->bic; 
	double new_aic = INFINITY, pre_aic = mod->aic;
	double JC_ll = 0;
	int false_positive;
	int diagnose = 0;	/* turn on diagnostic mode */
	unsigned int curr_K = K - 1;
	unsigned int n_param =  K * dat->max_read_length + curr_K // haplotypes and pi 
		+ (opt->use_error_profile?dat->n_quality*12:0);    // Previous is dat->n_quality * 16 

	mod->n_param = n_param;
	*fp = 1;

	if (!final) {

		/* check for indel errors */
		if (opt->nw_align == ALIGNMENT_HAPLOTYPES) {
			if ((err = check_fp_with_indels(opt, dat, mod, ini,
					curr_K, low_bound, error_profile,
							&false_positive)))
			 	return err;

			debug_msg(DEBUG_I, fxn_debug, "indel errors?: %i.\n",
								false_positive);

			if (false_positive)
				return NO_ERROR;  // *fp = 1
		}

		/* contamination diagnostic */		
		unsigned int rlen = ini->seed_lengths[curr_K];
		unsigned char *rseq = ini->seeds[curr_K];
		size_t *idx_array = NULL;
		unsigned int count = ini->uniq_seq_count[ini->H[curr_K]];

		if ((err = find_index(dat->seq_count, rseq, rlen, &idx_array)))
			return err;

		if ((err = abun_pvalue(opt, ini, idx_array, ini->e_trans,
				count, curr_K, opt->contamination_threshold, &ini->H_pvalue[curr_K],
							dat->sample_size, 0)))
			return err;
		
		/* control the type I error here */
		double adjp = opt->p_threshold;

		if (ini->H_pvalue[curr_K] > adjp) {

			debug_msg(DEBUG_I, fxn_debug, "remove false positive "
				"with Diagnostic Probability %8.2e with threshold %8.2e \n",
						ini->H_pvalue[curr_K], adjp);
			if (!diagnose)
				return NO_ERROR; // *fp = 1
		}

		if ((err = mean_exp_errors(dat, ini->uniq_seq_idx[ord],
				ini->uniq_seq_count[ord], &ini->H_ee[curr_K])))
			return err;

		debug_msg(DEBUG_II, fxn_debug, "The mean exp. n.err of %d is "
						"%.3f\n", ord, ini->H_ee[curr_K]);
	}

	/* recalculate transition probabilities here (currently only for
	 * ALIGNMENT_HAPLOTYPES) with new candidate haplotype
	 */
	/* [XP] I change to here from the main function haplotype_selection() because those transition
	 * [XP] probabilities are needed for bic check.  Since we may not get
	 * [XP] here, it can save time if it is a false positive by removed indel or
	 * [XP] contamination test.
	 */
	/* [XP] Other settings of options::nw_align are handled elsewhere;
	 * [XP] need more thought.
	 */
	/* Transition prob with alignment free strategy */
	if (opt->nw_align == ALIGNMENT_HAPLOTYPES && !final)
		ExpTrans_nogap(dat, opt, ini, ord, curr_K, error_profile,
							mod->err_encoding);

	/* update haplotypes */
	/* [TODO] modify the function and not use mod->haplotypes (just for Acceleration) */
	//for (unsigned int k = 0; k < K; k++)	/* [KSD] Haven't the first K-1 already been copied? */
	//memcpy(&mod->haplotypes[curr_K * dat->max_read_length],
	//		ini->seeds[curr_K], dat->max_read_length
	//					* sizeof *mod->haplotypes);

	/* update pi */
	Est_pi(ini, mod->pi, dat->sample_size, K, 1);

	/* calculate ll */
	mod->ll = Simple_Estep(mod, dat->sample_size, ini->e_trans, K);

	/* calculate aic and bic */
	if (opt->JC69_model) {
		modified_ic(ini->seeds[0], mod->est_ancestor, mod->distance,
				mod->ll, K, &JC_ll, &new_aic, &new_bic, n_param,
					dat->max_read_length, dat->sample_size,opt->ignor_nc);
	} else {
		new_aic = aic(mod->ll, n_param);
		new_bic = bic(mod->ll, n_param, dat->sample_size);
		JC_ll = 0;
	}

	mod->JC_ll = JC_ll;

//	if (final)
//		debug_msg(DEBUG_I, fxn_debug, "Final ll: %f; Final JC_ll: %f.\n",
//							mod->ll, mod->JC_ll);
//	else
		debug_msg(DEBUG_I, fxn_debug, "ll: %f; JC_ll: %f.\n",
							mod->ll, mod->JC_ll);

	/* better bic means the chosen seed is not a false positive */
	if (!diagnose && (opt->use_aic ? new_aic : new_bic) < (opt->use_aic ? pre_aic : pre_bic)) {
		*fp = 0;
		mod->bic = new_bic;
		mod->aic = new_aic;
	} else if (diagnose && (opt->use_aic ? new_aic : new_bic) < (opt->use_aic ? pre_aic : pre_bic)) {
		debug_msg(DEBUG_I, fxn_debug, "Would remove false positive by delta(BIC) = %f.\n", new_bic - pre_bic);
	}

	//if (!final) {
		debug_msg(DEBUG_I, fxn_debug, "pre_bic : %f. new_bic: %f (%f)\n",
					pre_bic, new_bic, new_bic - pre_bic);
		debug_msg(DEBUG_I, fxn_debug, "pre_aic : %f. new_aic: %f (%f)\n",
					pre_aic, new_aic, new_aic - pre_aic);
	//} else {
	//	debug_msg(DEBUG_I, fxn_debug, "Final bic: %f\n", new_bic);
	//	debug_msg(DEBUG_I, fxn_debug, "Final aic: %f\n", new_aic);
	//}
	

	return err;
}/* evaluate_haplotype */

/**
 * Estimate pi by:
 * (1) hard-clustering of reads assigned by maximum transition probability from
 *     each of K current haplotypes if \par reassign, or
 * (2) normalization of current estimated scaled true abundances if not
 *     \par reassign.
 * 
 * @param ini		pointer to initializer object
 * @param pi		cluster relative abundance (to be computed)
 * @param sample_size	total num of reads
 * @param K		current number of clusters
 * @param reassign	use hard-clustering by transition probability
 * 
 * @return		error status
 * 
 **/
int Est_pi(initializer *ini, double *pi, size_t sample_size,
					unsigned int K, int reassign)
{

	double max;
	unsigned int class;

	if (reassign) {
		for (unsigned int k = 0; k < K; k++)
			ini->cluster_size[k] = 0;

		for (unsigned int i = 0; i < sample_size; i++) {
			max = -INFINITY;
			class = 0;
			for (unsigned int k = 0; k < K; ++k) {
				if (max < ini->e_trans[k * sample_size + i]) {
					max = ini->e_trans[k * sample_size + i];
					class = k;
				}
			}
			ini->cluster_size[class]++;
			ini->cluster_id[i] = class;
		}

		for (unsigned int k = 0; k < K; ++k){
			pi[k] = (double) ini->cluster_size[k] / sample_size;
			pi[k] = log(pi[k]);
		}
	} else {
		double abun_sum = 0;
		for (unsigned int k = 0; k < K; k++)
			abun_sum += ini->H_abun[k];

		for (unsigned int k = 0; k < K; k++){
			pi[k] = ini->H_abun[k] / abun_sum;	
			pi[k] = log(pi[k]);
		}
	}

	return NO_ERROR;
}/* Est_pi */

/**
 * update seeds table in initializer when select a new haplotype
 * 
 * @param ini   			pointer to initializer object
 * @param dat   			pointer to data object
 * @param select            num of preselected haplotypes
 * @param ord 				idx of the new haplotype 
 * 
 * @return err    			error status
 * 
 **/
int update_seeds(data *dat, initializer *ini, unsigned int select, unsigned int ord){

	//ini->seed_idx[select] = ini->uniq_seq_idx[ord];  // idx in dmat and qmat
	ini->seed_lengths[select] = dat->lengths[ini->uniq_seq_idx[ord]];
	memcpy(ini->seeds[select], dat->dmat[ini->uniq_seq_idx[ord]],
		dat->max_read_length * sizeof **ini->seeds);
	ini->H[select] = ord ; // idx in unique sequence table
	ini->H_abun[select] = ini->abun_true[ord];   // estimated true abundance of haplotypes

	return NO_ERROR;
}/* update_seeds */


/**
 * Check candidate haplotype for possible indel misread of existing haplotypes.
 * Specifically, align the candidate haplotype to all previous haplotypes and
 * reestimate the scaled true abundance.  If it drops below the threshold,
 * the candidate haplotype is considered a "false positive".
 * 
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param low_bound	lower bound of haplotype abundance
 * @param error_profile	input error profile
 * @param select	number of haplotypes, excluding candidate
 * @param fp		detect as false positive (computed)
 * 
 * @return err		error status
 * 
 * * [TODO] it is too complex. Simplify it 
 **/
int check_fp_with_indels(options *opt, data *dat, model *mod, initializer *ini,
	unsigned int select, double low_bound, double *error_profile, int *fp)
{
	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;

	/* nw alignment for the candidate sequence and existing haplotypes */
	unsigned char ***nw_result = NULL;
	size_t *nw_alen = NULL;
	unsigned int *nw_mismatch = NULL;
	unsigned int *nw_indels = NULL;

	*fp = 0;

	nw_result = malloc(select * sizeof *nw_result);
	nw_alen = malloc(select * sizeof *nw_alen);
	nw_mismatch = malloc(select * sizeof *nw_mismatch);
	nw_indels = malloc(select * sizeof *nw_indels);

	if (!nw_result || !nw_alen || !nw_mismatch || !nw_indels) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "checkfp.nw");
		goto EXIT_CHECK_FP_WITH_INDELS;
	}

	int ends_free = 1;  // ends-free alignment: no penalty for terminal gaps

	/* haplotype candidate */
	unsigned int rlen = ini->seed_lengths[select];
	unsigned char *rseq = ini->seeds[select];

	/* align candidate haplotype to each existing haplotype */
	/* [KSD] duplicative code with nwalign_matrix() */
	for (unsigned int k = 0; k < select; k++) {
		size_t alen;
		unsigned int ref_len = ini->seed_lengths[k];
		unsigned char *ref_seq = ini->seeds[k];

		unsigned char **aln = nwalign(ref_seq, rseq, (size_t) ref_len,
			(size_t) rlen, opt->score, opt->gap_p, opt->band,
						ends_free, NULL, &err, &alen);


		#if DEBUG
		//fprintf(stderr, "alignment");
		for (size_t j = 0; j < alen; ++j) {
				fprintf(stderr, "%c", aln[0][j] == '-'
				? '-' : xy_to_char[(int) aln[0][j]]);
			}
			fprintf(stderr, "\n");
			for (size_t j = 0; j < alen; ++j) {
				fprintf(stderr, "%c", aln[1][j] == '-'
				? '-' : xy_to_char[(int) aln[1][j]]);
			}
			fprintf(stderr, "\n");
		#endif

		nw_result[k] = aln;
		nw_alen[k] = alen;

		/* calculate number of indels and mismatch based on alignment */
		ana_alignment(aln, alen, rlen, &nw_indels[k], 
					&nw_mismatch[k], opt->ends_free,opt->info);

		/* diagnostic output: */
		//debug_msg(DEBUG_I, DEBUG_I, "haplotype %u alignment length %zu; #indels %u; #mismatches %u\n", k, alen, nw_indels[k], nw_mismatch[k]);
	}

	/*------------------------------------------------------------------- */
	/* recalculate transition prob for each read with same sequence as
	 * candidate haplotype now using NW alignment; need to allocate memory
	 */

	double *e_trans = NULL;
	size_t *idx_array = NULL;
	unsigned int count = ini->uniq_seq_count[ini->H[select]];

	if ((err = find_index(dat->seq_count, rseq, rlen, &idx_array)))
		goto EXIT_CHECK_FP_WITH_INDELS;

	e_trans = malloc(select * count * sizeof *e_trans);

	if (!e_trans) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "checkfp.e_trans");
		goto EXIT_CHECK_FP_WITH_INDELS;
	}

	for (unsigned int k = 0; k < select; k++) {
		double tmp = 0;

		for (unsigned int r = 0; r < count; r++) {
			e_trans[k*count + r] = trans_nw(opt, nw_result[k],
				nw_alen[k], nw_mismatch[k], nw_indels[k],
				error_profile, mod->err_encoding,
				dat->qmat[idx_array[r]], dat->n_quality, 
				mod->adj_trunpois, rlen, dat->error_prob,opt->ends_free);
			tmp += e_trans[k*count + r];
			#if DEBUG
			if (nw_indels[k] == 1) {
				//mmessage(INFO_MSG, NO_ERROR, "indel rate: %f\n", temp);
				mmessage(INFO_MSG, NO_ERROR, "the trans prob of the %i th read to hap %i is %f \n",r, k, e_trans[k*count + r]);
				mmessage(INFO_MSG, NO_ERROR, "The self trans is %f\n", ini->self_trans[idx_array[r]]);
			}
			#endif
		}
	}

	/* ------------------------------------------------------------------ */
	/* reestimate scaled true abundance of candidate haplotype            */
	/* [KSD] duplicative code with other locations; merge to avoid
	 * [KSD] introducing differences in, e.g., convergence criteria.
	 */

	double true_abun = count;
	double sum_pro_H, true_abun_new;
	unsigned int iter = 0;
	double delta;

	do {
		true_abun_new = count;
		//true_abun_var = 0;

		for (unsigned int r = 0; r < count; r++) {

			double self_ts = ini->self_trans[idx_array[r]];

			/* modify self trans rate if we consider indels */
			if (opt->indel_model == INDEL_PER_READ)
				self_ts += dpois(0, opt->indel_error * rlen, 1) - mod->adj_trunpois;

			sum_pro_H = 0.;
			for (unsigned int k =0 ; k < select; k++)
				sum_pro_H += exp(e_trans[k * count + r])
							* ini->H_abun[k];
			true_abun_new -= sum_pro_H
				/ (sum_pro_H + exp(self_ts) * true_abun);
		}

		delta = (true_abun - true_abun_new) / true_abun;

		if (delta < 0 || true_abun_new < 0)
			return mmessage(WARNING_MSG, INTERNAL_ERROR,
				"True abundnace increase or under 0\n");

		if (delta < opt->epsilon || true_abun_new < low_bound)
			break;

		true_abun = true_abun_new;
		iter ++;	

		if (iter > opt->n_iter_amplici)
			return mmessage(WARNING_MSG, INTERNAL_ERROR,
				"exceed max interations,%i\n", iter);

	} while (1);

	/* just a simple check for indel errors */
	if (true_abun_new < low_bound)
		*fp = 1;

	debug_msg(DEBUG_I, fxn_debug,"observe abundance: %i; previous "
				"abundance: %f; new_abundance: %f \n", count,
					ini->H_abun[select], true_abun_new);
	if (!(*fp))
		/* use true_abun_new as the current abundance */
		ini->H_abun[select] = true_abun_new;

EXIT_CHECK_FP_WITH_INDELS:

	if (nw_alen)
		free(nw_alen);
	if (nw_result){
		free_nw_result(nw_result,select);
		free(nw_result);
	}
	if (nw_mismatch)
		free(nw_mismatch);
	if (nw_indels)
		free(nw_indels);
	if (e_trans)
		free(e_trans);

	return err;
}/* check_fp_with_indels */

/**
 * Contamination diagnosis.  We perform a simple disgnostic test for
 * contamination, where the contamination process is assumed to produce a fixed
 * number of reads of every candidate.  If the number of observed reads is
 * considerably higher than this threshold and not explained by misreads from
 * the current haplotype set, then we believe the candidate haplotype is real.
 * 
 * @param opt		pointer to options object
 * @param ini		pointer to initializer object
 * @param idx_array	indices of reads matching candidate haplotype
 * @param count		abserved abundance of candidate haplotype
 * @param e_trans	transition probability matrix
 * @param select	index of candidate haplotype
 * @param threshold  threshold for contamination
 * @param error_profile	input error profile
 * @param p		pointer to p-value
 * @param sample_size	total number of reads
 * @param partial	is e_trans for all data or just the reads matching
 *			candidate haplotype
 * 
 * @return err 		error status
 * 
 **/
int abun_pvalue(options *opt, initializer *ini, size_t *idx_array,
	double *e_trans, unsigned int count, unsigned int select,
	unsigned int threshold, double *p, size_t sample_size, int partial)
{

	int fxn_debug = opt->info;

	double true_abun_var = 0.;
	double abun_null = 0.;  //  estimated number of count from others under null hypothesis
	double gamma = 0.;
	unsigned int rlen = ini->seed_lengths[select];
	double true_abun = ini->H_abun[select];  // note the abundance of the haplotype is fixed when generate p value
	double *perr = NULL;

	/* compute exact p value if n is small enough */
	if (count < MAX_N_EXACT_P) {
		perr = malloc(count * sizeof *perr);
		if (!perr)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "perr");
	}

	for (unsigned int r = 0; r < count; r++) {

		double self_ts = ini->self_trans[idx_array[r]];

		/* modify self trans rate if we consider indels (no longer use) */ 
		//if (opt->indel_model == INDEL_PER_READ)
		//	self_ts += dpois(0, opt->indel_error * rlen, 1);

		double Ey = exp(self_ts) * true_abun;   //previous  double Ey = exp(self_ts)*true_abun
		double Ex = 0.;
		for (unsigned int k =0 ; k < select; k++) {
			if (partial)
				Ex += exp(e_trans[k * count + r]) * ini->H_abun[k];
			else
				Ex += exp(e_trans[k * sample_size + idx_array[r]])
							* ini->H_abun[k];
		}
		//mmessage(INFO_MSG,NO_ERROR,"Ey:%f; Ex: %f \n", Ey, Ex);
		double p = Ex / (Ex + Ey);
		double Varf = p * (1 - p);
		if (perr)
			perr[r] = p;
		//mmessage(INFO_MSG,NO_ERROR,"%f  \n",p);

		true_abun_var += Varf;
		abun_null += p;
		gamma += p * (1 - p) * (1 - 2 * p);
	}

	debug_msg(DEBUG_II, fxn_debug, "variance under null: %f \n", true_abun_var);
	debug_msg(DEBUG_II, fxn_debug, "N_H->sm under null: %f \n", abun_null);
	
	int bound = count - (threshold+1); 
	debug_msg(DEBUG_I, fxn_debug, "bound=%i;\n", bound);
	debug_msg(DEBUG_I, fxn_debug, "count=%i;\n", count);

	if (perr){
		*p = ppoisbin(bound, count, perr, 1); // P(S > bound)
	}else{ /* approximate p value */
		// *p = ppois(lower_bound, true_abun, 1, 0);	/* [KSD] What is this? */
		/* avoid numeric problem here */
		double sigma = sqrt(true_abun_var);
		/***
		double x;
		if (fabs(sigma-0.) < DBL_EPSILON) {
			*p = 0;
			goto EXITNOW;
			gamma = INFINITY;
			x = INFINITY;
		} else {
			gamma = gamma / (sigma * sigma * sigma);
			x = (bound + 0.5 - abun_null) / sigma;
		}
		double phi_x = dnorm(x,0,1,0);
		//double cdf_x = pnorm(x,0,1,0,0);
		debug_msg(DEBUG_II, fxn_debug, "gamma: %f,x: %f, phi_x:%f \n",
							gamma, x, phi_x);
		if (fabs(phi_x - 0.) < DBL_EPSILON)
			*p = pnorm(x,0,1,1,0);
		else
			*p = pnorm(x,0,1,1,0) + phi_x*gamma*(1.0-x*x)/6; 
		*p = 1.0 - *p; // 1- F(x) 
		if(*p < 0) // avoid bug
			*p = 0.;
		 //refined normal approximation 
		 ***/
		 *p = pnorm(bound, abun_null, sigma, 0, 0); // normal approximation
		//*p = ppois(bound,abun_null,0,0); //poisson approximation
	}
	debug_msg(DEBUG_I, fxn_debug, "Diagnostic Probability=%8.2e;\n", *p);

	if (*p < 0 || *p >1)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
				"Diagnostic Probability is not in [0,1] ");
EXITNOW:

	if (perr)
		free(perr);

	return NO_ERROR;
} /* abun_pvalue */

/**
 * Compute model log likelihood given estimates of pi and computed
 * transition probabilities from each existing haplotype.  Also
 * computes posterior probabilities of assignment for each read
 * using a underflow-safe calculation.
 * 
 * @param mod		pointer to model object
 * @param sample_size	total number of reads
 * @param e_trans	transition probability matrix	
 * @param K		number of haplotypes
 * 
 * @return ll		model log likelihood
 **/
double Simple_Estep(model *mod, size_t sample_size, double *e_trans,
							unsigned int K)
{
	
	double ll = 0;
	double max, sum;
	unsigned int idx;

	for (unsigned int i = 0; i < sample_size; i++) {
		max = -INFINITY;
		for (unsigned int k = 0; k < K; k++) {
			idx = k * sample_size + i;
			mod->eik[idx] = mod->pi[k] + e_trans[idx];
			if (max < mod->eik[idx])
				max = mod->eik[idx];
		}
		sum = 0.;
		for (unsigned int k = 0; k < K; ++k) {
			idx = k * sample_size + i;
			mod->eik[idx] = exp(mod->eik[idx] - max);
			sum += mod->eik[idx];
		}

		ll += log(sum) + max;

		/* actually these posterior probabilities are not used */
		for (unsigned int k = 0; k < K; ++k)
			mod->eik[k * sample_size + i] /= sum;
			
	}

	return ll;
}/* Simple_Estep */

/* print reads assignment */
void fprint_assignment(FILE *fp, unsigned int *v, size_t n, unsigned int max, int width, int newline){
	size_t i;
	for (i = 0; i < n; ++i){
		if(v[i]< max){
			if (width)
				fprintf(fp, " %*u", width, v[i]);
			else
				fprintf(fp, " %u", v[i]);
		}else
			fprintf(fp, " NA");
			
	}
	if (newline) fprintf(fp, "\n");
} /* fprint_assignment */

/**
 * Expected number of errors for one read.
 * 
 * @param qual		quality score seq
 * @param length	length of the seq
 * @param error_prob	error prob indicated by quality score
 * 
 **/
double exp_errors(unsigned char *qual, unsigned int length, double *error_prob)
{
		
	double E_err = 0.;

	for (unsigned int j = 0; j < length; j++)
		E_err += error_prob[qual[j]];

	return E_err;
}/* exp_errors */

/**
 * Compute mean expected number of errors for each uniq seq.  [KSD] Why?
 * 
 * @param dat		pointer to data object
 * @param idx		idx of the uniq seq in dmat
 * @param count_i	observed abundance of the uniq seq
 * @param mean_exp_err	pointer to mean_exp_err to set
 * @return		error status
 **/
int mean_exp_errors(data *dat, unsigned int idx, unsigned int count_i,
							double *mean_exp_err)
{
	int err;

	unsigned char *seq = dat->dmat[idx];
	unsigned int length = dat->lengths[idx];

	size_t *idx_array = NULL;
	if (((err = find_index(dat->seq_count, seq, length, &idx_array))))
		return err;

	double temp;
	double sum = 0.;

	for (unsigned int i = 0; i < count_i; i++) {
		temp = exp_errors(dat->qmat[idx_array[i]], length,
							dat->error_prob);
		sum += temp;
	}

	*mean_exp_err = sum / count_i;

	return err;
}/* mean_exp_errors */


/* output the format fasta file for UCHIME (use size) */
void fprint_haplotypes_size(FILE *fp, data_t *data, size_t n, size_t p, double pthres, 
			char const * const prefix, double *pvalue, unsigned int *size, double *ee) {
	for (size_t i = 0; i < n; ++i) {
		if(pvalue && pvalue[i]>=pthres)
			continue;
		fprintf(fp, ">%s%lu;", prefix, i);
		if(size) fprintf(fp,"size=%u;",size[i]);
		if(pvalue) fprintf(fp,"DiagP=%8.2e;",pvalue[i]);
		if(ee) fprintf(fp, "ee=%.3f;",ee[i]);
		fprintf(fp, "\n");
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i*p + j]]);
		fprintf(fp, "\n");
		}
} /* fprint_haplotypes_size */

/* output the format fasta file for UCHIME (use relative true abundance */
void fprint_haplotypes_abun(FILE *fp, data_t *data, size_t n, size_t p, double pthres, char const * const prefix, 
							double *pvalue, double *abun, double *ee) {
	for (size_t i = 0; i < n; ++i) {
		if(pvalue && pvalue[i]>=pthres)
			continue;
		fprintf(fp, ">%s%lu;", prefix, i);
		if(abun) fprintf(fp,"size=%.3f;",abun[i]);
		if(pvalue) fprintf(fp,"DiagP=%8.2e;",pvalue[i]);
		if(ee) fprintf(fp, "ee=%.3f;",ee[i]);
		fprintf(fp, "\n");
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i*p + j]]);
		fprintf(fp, "\n");
	}
} /* fprint_haplotypes_abun */

/**
 * Assign reads to clusters while filtering on log likelihood or posterior
 * probability.  If a read is filtered out, it gets assigned to cluster with
 * index K, aka NA.
 * 
 * @param K		number of clusters
 * @param ll_cutoff	cutoff of log likelihood to assign read
 * @param eik		maximum posterior assignment probability (log)
 *			if eik ignore pi and e_trans
 * @param pi		relative abundance for each clusters (log)
 * @param e_trans	transition prob matrix
 * @param sample_size	total number fo reads
 * @param ri		pointer to run_info object
 * 
 * @return err		err status
 **/
int likelihood_filter(unsigned int K, double ll_cutoff, double *eik, double *pi,
	double *e_trans, size_t sample_size, run_info *ri)
{
	unsigned int bound = K + 1; 

	for (unsigned int k = 0; k < bound; k++)
		ri->optimal_cluster_size[k] = 0;

	/* compute & store maximum conditional log likelihood */
	for (unsigned int i = 0 ; i < sample_size; i++) {
		if (eik)
			ri->optimal_cluster_ll[i] = eik[ri->optimal_cluster_id[i] * sample_size + i];
		else
			ri->optimal_cluster_ll[i] = pi[ri->optimal_cluster_id[i]]
				+ e_trans[ri->optimal_cluster_id[i] * sample_size + i];

		//ri->optimal_cluster_ll[i] = eik[ri->optimal_cluster_id[i] * sample_size + i];  

		if (ri->optimal_cluster_ll[i] < ll_cutoff)
			ri->optimal_cluster_id[i] = K;         // treat those outliers as a new cluster

		ri->optimal_cluster_size[ri->optimal_cluster_id[i]]++;
	}

	return NO_ERROR;
}/* likelihood_filter */

/**
 * get MMEs of all parameters in the JC69 model. 
 * 
 * @param hap	pointer to haplotypes
 * @param anc	ancestor haplotype, to be calculated (tbc)
 * @param dist	expected no. changes/site b/w ancestor & haplotype, tbc
 * @param K	number of haplotypes
 * @param len	length of reads
 * @return err
 **/
int m_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	unsigned int K, unsigned int len, int start)
{
	
	int err = NO_ERROR;
	unsigned int count[NUM_NUCLEOTIDES];
	unsigned int max_count;
	// int start = 9;  // ignore the first 9 nucleotides

	/* most common nucleotide across haplotypes is the estimated ancestor */
	for (unsigned int j = start; j < len; j ++){
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; n++)
			count[n] = 0;
		for (unsigned int k = 0; k < K; k++)
			count[hap[k * len + j]]++;
		max_count = 0;
		for (unsigned char n = 0; n < NUM_NUCLEOTIDES; ++n){
			if (count[n] > max_count) {
				max_count = count[n];
				anc[j] = n; 
			}
		}
	}
	for (unsigned int j = 0; j < start; j ++)
		anc[j] = (unsigned char) 4;


	/* estimate the expected number of changes per site of all haplotypes */
	for (unsigned int k = 0; k < K; k++) {		
		double tmp = (double) hamming_char_dis( (char *) &hap[k*len+start],
					(char *) anc, (size_t) len-start) / (len-start);
		/* Previous there is a bug for the estimated distance out of range */
		if(tmp >= 0.75)
			dist[k] = INFINITY;
		else
			dist[k] = -0.75 * log(1 - tmp / 0.75);
	}

	return err;
}/* m_JC69 */


/**
 * Calculate the log likelihood of all K haplotypes under JC69 model.
 * 
 * @param hap   haplotype sequences
 * @param anc   ancestor sequence
 * @param dis	expected no. changes/site for each haplotype
 * @param K	number of haplotypes
 * @param len	length of the haplotypes
 * @return	err status
 **/
double e_JC69(unsigned char * hap, unsigned char * anc, double *dist,
	unsigned int K, unsigned int len,int start)
{
	double ll = 0;
	// int start = 9;  // ignore the first 9 nucleotides

	for (unsigned int k = 0; k < K; k++)
		for (unsigned int j = start; j < len; j++)
			if (anc[j] == hap[k * len + j])
				ll += log(0.25 + 0.75 * exp(-dist[k] / 0.75));
			else
				ll += log(0.25 - 0.25 * exp(-dist[k] / 0.75));

	return ll;
}/* e_JC69 */

/**
 * Calculate aic and bic modified by approximate JC69 hierarchical model on
 * haplotypes.
 *
 * @param hap			haplotypes
 * @param est_anc		ancester sequence
 * @param distance		distance from haplotypes to ancestor sequence
 * @param best_ll		current log likelihood from data
 * @param K			number of haplotypes
 * @param JC_ll			log likelihood from JC69 model
 * @param n_aic			pointer to aic, value updated
 * @param n_bic			pointer to bic, value updated
 * @param n_param		number of parameters in current model
 * @param max_read_length	length of haplotypes
 * @param sample_size		sample size
 *
 * return			err status
 **/
int modified_ic(unsigned char *hap, unsigned char *est_anc, double *distance,
	double best_ll, unsigned int K, double *JC_ll, double *n_aic,
	double *n_bic, unsigned int n_param, unsigned int max_read_length,
	size_t sample_size, int start)
{
	int param_change = 0;
	// int start = 9;  // ignore the first 9 nucleotides

	
	m_JC69(hap, est_anc, distance, K, max_read_length,start);
	*JC_ll = e_JC69(hap, est_anc, distance, K, max_read_length,start);

	/* K branch lengths, ancestral haplotype, but no haplotypes estimated */
	 param_change = K - (max_read_length -start) * (K - 1);

	*n_aic = aic(best_ll + *JC_ll, n_param + param_change);
	*n_bic = bic(best_ll + *JC_ll, n_param + param_change, sample_size);

	return NO_ERROR;
}/* modified_ic */

/**
 * Reads assignment with given haplotype set 
 * 
 * @param ini	pointer to initializer object
 * @param dat	pointer data object
 * @param opt	pointer to options object
 * @param mod   pointer to model object
 * @param ri    pointer to run_info object
 *
 * return	err status
 **/
int reads_assignment(options * opt, data * dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;
	int fxn_debug = opt->info;
	double l1third = 1./3;

	/* maybe use error profile */
	double *error_profile = NULL;
	if (opt->use_error_profile && mod->error_profile) {
		error_profile = mod->error_profile;
		debug_msg(DEBUG_II, fxn_debug, "Use error profile. \n");
	}

	if((err = trans_expectation(opt, dat, ini, error_profile, 
					mod->adj_trunpois, mod->eik,0)))
		return err;
	/* Keep codes below for debug purpose  */
	/*
	
	for(unsigned int u = 0; u <dat->hash_length; ++u ){
		unsigned char *read = dat->dmat[ini->uniq_seq_idx[u]];
		unsigned int rlen = dat->lengths[ini->uniq_seq_idx[u]];

		unsigned int count = ini->uniq_seq_count[u]; // num. of reads wih unique seq
		size_t *idx_array;  // idx of reads

		if ((err = find_index(dat->seq_count, read, rlen, &idx_array)))
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Cannot find in the hash table !");

		//  align to haplotypes  
		for (unsigned int h = 0; h < opt->K; ++h){

			unsigned char *hap_seq = ini->seeds[h];

			if(opt->nw_align == NO_ALIGNMENT){
				
				for(unsigned int r = 0; r < count; ++r){
					double eik = 0.;
					size_t id = idx_array[r];
					for (unsigned int j = 0; j < dat->lengths[r]; j++) {

						if (error_profile) {
							if (mod->err_encoding == STD_ENCODING)
								eik += translate_error_STD_to_XY(
									error_profile,
									dat->n_quality, hap_seq[j],
									dat->dmat[id][j],
									dat->qmat[id][j]);
						else if (mod->err_encoding == XY_ENCODING)
							eik += error_profile[(NUM_NUCLEOTIDES
								* hap_seq[j] + dat->dmat[id][j])
								* dat->n_quality
								+ dat->qmat[id][j]];
						} else {
							//double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);
							double ep = dat->error_prob[dat->qmat[id][j]];	
							if (dat->dmat[id][j] == hap_seq[j] )			
								eik += log(1 - ep);
							else
								eik += log(ep) + l1third;
						}
					}
					mod->eik[h*dat->sample_size+ id] = eik;		
				}

			}else{
				size_t alen = dat->max_read_length;
				unsigned int nindels = 0;
				unsigned int nmismatch = 0;

				unsigned char **aln = nwalign(hap_seq, read,
				(size_t) ini->seed_lengths[h],
				(size_t) rlen,
				opt->score, opt->gap_p, opt->band, 1, NULL,
									&err, &alen);

				// count for number of indels 
				ana_alignment(aln, alen, rlen, &nindels, 
						&nmismatch, opt->info);

				for(unsigned int r = 0; r<count;++r){

					mod->eik[h*dat->sample_size+ idx_array[r]] = trans_nw(opt, aln,
						alen, nmismatch, nindels, error_profile, mod->err_encoding,
						dat->qmat[idx_array[r]], dat->n_quality, mod->adj_trunpois,
										rlen, dat->error_prob);
					
					debug_msg(DEBUG_III, fxn_debug, "num of indels: %i; num of "
							"mismatch: %i\n", nindels, nmismatch);

				}
				// free 
				if(aln){
					free(aln[0]);
					free(aln[1]);
					free(aln);
					aln = NULL;
				}
			}
		}
	}
	*/

	if(opt->trans_matrix){
		FILE *fp = fopen(opt->trans_matrix, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->trans_matrix);
		//fprint_vectorized_matrix(fp,mod->eik,dat->sample_size,opt->K,0);  not work 
		for (size_t i = 0; i < dat->sample_size; ++i) {
		    fprintf(fp, "%3lu", i);
			for (unsigned int j = 0; j < opt->K; ++j) {
				fprintf(fp, " %8.2e", mod->eik[j*dat->sample_size + i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	/* simply update mod->pi */
	assign_clusters(mod->eik, opt->K, dat->sample_size, ri->optimal_cluster_size,
		ri->optimal_cluster_id, 1);
	for (unsigned int k = 0; k < opt->K; ++k) {
		mod->pi[k] = (double) ri->optimal_cluster_size[k] / dat->sample_size;
		if (!mod->pi[k])
			mod->pi[k] = 1.0 / dat->sample_size;  // possible if given haplotypes 
		mod->pi[k] = log(mod->pi[k]);
	}
	
	/* update mod->eik with new estimated mod->pi */
	for (unsigned int r = 0; r<dat->sample_size; ++r)
		for (unsigned int k = 0; k < opt->K; ++k)
			mod->eik[k*dat->sample_size+r] += mod->pi[k];

	/* reassign reads with updated mod->eik (unnormalized ) */
	assign_clusters(mod->eik, opt->K, dat->sample_size, ri->optimal_cluster_size,
		ri->optimal_cluster_id, 1);

	/* filter with unnormalized mod->eik (pi* e_trans) */
	likelihood_filter(opt->K, opt->ll_cutoff, mod->eik, NULL, NULL,
							dat->sample_size, ri);

	/* output the reads assignment */
	FILE *fp = NULL;

	opt->outfile_info = opt->outfile_base;

	fp = fopen(opt->outfile_info, "w");
    if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->outfile_info);

	fprintf(fp, "assignments: ");
	fprint_assignment(fp, ri->optimal_cluster_id, dat->sample_size,
								opt->K, 2, 1);
	fprintf(fp, "cluster sizes: ");
	fprint_uints(fp, ri->optimal_cluster_size, opt->K, 3, 1);

	fclose(fp);
		
	mmessage(INFO_MSG, NO_ERROR, "Output the assignment"
				"file: %s \n", opt->outfile_info);

	return err;
} /* reads_assignment */

/* calculate transition probability between reads and haplotypes */
int trans_expectation(options *opt, data *dat,initializer*ini, double *error_profile, 
					double adj_trunpois, double *trans_prob, int ends_free){

	int err = NO_ERROR;
	double l1third = 1./3;
	int fxn_debug = opt->info;

	for(unsigned int u = 0; u < dat->hash_length; ++u ){

		//mmessage(INFO_MSG, NO_ERROR, "The %5d uniq seq of total %5d sequences\n", u, dat->hash_length);

		unsigned char *read = dat->dmat[ini->uniq_seq_idx[u]];
		unsigned int rlen = dat->lengths[ini->uniq_seq_idx[u]];

		unsigned int count = ini->uniq_seq_count[u]; // num. of reads wih unique seq
		size_t *idx_array;  // idx of reads

		if ((err = find_index(dat->seq_count, read, rlen, &idx_array)))
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Cannot find in the hash table !");

		/* align to haplotypes  */
		for (unsigned int h = 0; h < opt->K; ++h){

			unsigned char *hap_seq = ini->seeds[h];

			if(opt->nw_align == NO_ALIGNMENT){
				
				for(unsigned int r = 0; r < count; ++r){
					double eik = 0.;
					size_t id = idx_array[r];
					for (unsigned int j = 0; j < dat->lengths[r]; j++) {

						if (error_profile) {
							if (opt->err_encoding == STD_ENCODING)  // make sure opt->err_encoding is the same as mod->err_encoding
								eik += translate_error_STD_to_XY(
									error_profile,
									dat->n_quality, hap_seq[j],
									dat->dmat[id][j],
									dat->qmat[id][j]);
						else if (opt->err_encoding == XY_ENCODING)
							eik += error_profile[(NUM_NUCLEOTIDES
								* hap_seq[j] + dat->dmat[id][j])
								* dat->n_quality
								+ dat->qmat[id][j]];
						} else {
							//double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);
							double ep = dat->error_prob[dat->qmat[id][j]];	
							if (dat->dmat[id][j] == hap_seq[j] )			
								eik += log(1 - ep);
							else
								eik += log(ep) + l1third;
						}
					}
					trans_prob[h*dat->sample_size+ id] = eik;		
				}

			}else{
				size_t alen = dat->max_read_length;
				unsigned int nindels = 0;
				unsigned int nmismatch = 0;

				unsigned char **aln = nwalign(hap_seq, read,
				(size_t) ini->seed_lengths[h],
				(size_t) rlen,
				opt->score, opt->gap_p, opt->band, 1, NULL,
									&err, &alen);

				/* count for number of indels */
				ana_alignment(aln, alen, rlen, &nindels, 
						&nmismatch, opt->ends_free,opt->info); // need further check 

				for(unsigned int r = 0; r<count;++r){

					/*	if(idx_array[r] ==0 && h == 2 ){
					for (size_t j = 0; j < alen; ++j) {
					fprintf(stderr, "%c", aln[0][j] == '-'
					? '-' : xy_to_char[(int) aln[0][j]]);
					}
					fprintf(stderr, "\n");
					for (size_t j = 0; j < alen; ++j) {
					fprintf(stderr, "%c", aln[1][j] == '-'
					? '-' : xy_to_char[(int) aln[1][j]]);
					}
					fprintf(stderr, "\n");
					}
					*/

					trans_prob[h*dat->sample_size+ idx_array[r]] = trans_nw(opt, aln,
						alen, nmismatch, nindels, error_profile, opt->err_encoding,
						dat->qmat[idx_array[r]], dat->n_quality, adj_trunpois,
										rlen, dat->error_prob,ends_free);
					/*
					if(idx_array[r] ==0 && h == 2 ){
						 for (size_t j = 0; j < alen; ++j){
                        	fprintf(stderr, " %8.2e\n", error_profile[(
							NUM_NUCLEOTIDES
							* aln[0][j] + aln[1][j])
							* dat->n_quality +dat->qmat[idx_array[r]][j]]);
							fprintf(stderr, " %c\n", dat->qmat[idx_array[r]][j]+dat->fdata->min_quality);
                    	}
					 fprintf(stderr, "\n");	
					fprintf(stderr, " %8.2e\n", trans_prob[h*dat->sample_size+ idx_array[r]]);
					}
					*/
					
					//debug_msg(DEBUG_III, fxn_debug, "num of indels: %i; num of "
					//		"mismatch: %i\n", nindels, nmismatch);

				}
				/* free */
				if(aln){
					free(aln[0]);
					free(aln[1]);
					free(aln);
					aln = NULL;
				}
				if(err)
					return err;
			}
		}
	}

	return err;
}