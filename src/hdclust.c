/**
 * @file hdclust.c
 * @author Karin S. Dorman
 *
 * Cluster greedily by Hamming distance, similar to UMI-tools.
 *
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <math.h>

#include "hdclust.h"
#include "ampliclust.h"
#include "initialize.h"
#include "statistics.h"
#include "io.h"
#include "hash.h"
#include "align.h"
#include "error.h"
#include "error_est.h"

/* amplici */
int hdclust_malloc(options *opt, data *dat, initializer *ini, unsigned int K_space, unsigned n_candidate);
int hdclust_realloc(options *opt, initializer *ini, model *mod, unsigned int preK, unsigned int K, size_t sample_size, unsigned int hash_length, unsigned int max_read_length);

/* estimated scaled true abundance */
int expected_TrueAbundance(options *opt, data *dat, double *H_abun, double *e_trans, double *self_trans, double *abun_true, size_t *idx, unsigned int count_i, unsigned int select, unsigned int i, int conve, double low_bound);
double iterate_expected_true_abundance(unsigned int sample_size, size_t *idx_array, double *e_trans, double *self_trans, double *H_abun, unsigned int select, unsigned int obs_abun, double true_abun);

/* check false positives */
int check_indel_error(options *opt, data *dat, model *mod, initializer *ini, unsigned int K, unsigned int ord);
int check_fp_with_indels(options *opt, data *dat, model *mod, initializer *ini, unsigned int select, double low_bound, double * error_profile, int *fp);
int Est_pi(initializer *ini, double *pi, size_t sample_size, unsigned int K, int reassign);
int abun_pvalue(options *opt, initializer *ini, size_t *idx_array, double *e_trans, unsigned int count, unsigned int select, unsigned int threshold, double *p, size_t sample_size, int partial);


/* transition prob with or without alignment */
int Expected_SelfTrans(options *opt, model *mod, data *dat, double *self_trans, double *error_profile, int err_encoding);
int ExpTrans_nogap(data *dat, options *opt, initializer *ini, unsigned int uidx, unsigned int select, double *error_profile, int err_encoding);
int ExpTrans_nwalign(model *mod, data *dat, options *opt, initializer *ini, data_t ***nw_result, size_t *nw_alen, unsigned int select, double *error_profile, int err_encoding, unsigned int H_id);

void fprint_haplotype(FILE *fp, data_t *data, unsigned int len);

/**
 * Cluster reads.
 *
 * @param ini	pointer to initializer object
 * @param dat	pointer data object
 * @param opt	pointer to options object
 * @param mod	pointer to model object
 * @param ri	pointer to run_info object
 *
 * return	err status
 **/
int hd_clust(options *opt, data *dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;
	int output_hap = 1;	/* output haplotype FASTA file */
	int use_size = 0;	/* output cluster sizes for abundance */

	/* K to be determined in the function below, opt->K will be changed */
	if ((err = cluster_by_hd(opt, dat, mod, ini, opt->hamming_proportion)))
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

		fprintf(fp, "pi: ");
		for(unsigned int k = 0; k < opt->K; ++k)
			mod->pi[k] = exp(mod->pi[k]);
		fprint_doubles(fp, mod->pi, opt->K ,6,1);

		fprintf(fp, "reads ll: ");
		fprint_doubles(fp, ri->optimal_cluster_ll, dat->sample_size,
									3, 1);

		//[TODO] output ini->seeds directly
		fprint_fasta(fp, ini->seeds[0], opt->K,
					 dat->max_read_length, ini->seed_lengths, "H");

		fprintf(fp, "ee: "); // mean expected number of errors
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
		if (opt->run_diagnostic_test) {
			fprintf(fp, "p value: ");
			for (unsigned k = 0; k < opt->K; k++){
				double pvalue = exp(ini->H_pvalue[k]);
				if (pvalue < 1e-3)
					fprintf(fp, " %8.2e", pvalue);
				else
					fprintf(fp, " %.3f", pvalue);
			}
			fprintf(fp, "\n");
		}
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

		/* two function allow variable length */
		if (use_size)
			fprint_haplotypes_size(fp2, ini->seeds, opt->K,
				ini->seed_lengths, opt->p_threshold, "H",
				opt->run_diagnostic_test ? ini->H_pvalue : NULL,
					ri->optimal_cluster_size, ini->H_ee);
		else
			fprint_haplotypes_abun(fp2, ini->seeds, opt->K,
				ini->seed_lengths, opt->p_threshold, "H",
				opt->run_diagnostic_test ? ini->H_pvalue : NULL,
							ini->H_abun, ini->H_ee);

		fclose(fp2);

		mmessage(INFO_MSG, NO_ERROR, "Output the final haplotype fasta "
					"file: %s \n", opt->outfile_fasta);

	}

	if (outfile_hap)
		free(outfile_hap);

	return err;
}/* hd_clust */

/**
 * Select no more than K_max real haplotype sequences
 *
 * @param ini	initializer object
 * @param dat	data object
 * @param opt	options object
 * @param mod	model object
 * @param hprop	maximum Hamming proportion
 *
 * return	err status
 **/
int cluster_by_hd(options *opt, data *dat, model *mod, initializer *ini,
								double hprop)
{
	int err = NO_ERROR;
#ifdef DEBUG_AMPLICI
	int fxn_debug = opt->info;//DEBUG_I;//
#endif

#ifdef DEBUG_AMPLICI
	debug_msg(DEBUG_I, fxn_debug, "hp=%f nu=%u q=[%u,%u]\n", hprop,
			dat->hash_length, dat->min_quality, dat->max_quality);
#endif

	/* ------------------------------------------------------------------ */
	/* Variable Declaration */

	unsigned int K_curr = 0;	/* current number of haplotypes */
	unsigned int curr_useq_idx = 0;	/* unique sequence index of current haplotype */

	unsigned int K_space = opt->K;	/* current space for clusters */
	unsigned int ini_K = opt->K;	/* the initial K space */
	unsigned int n_candidate = dat->hash_length;

#ifdef DEBUG_AMPLICI
	debug_msg(DEBUG_II, fxn_debug, "Number of candidates: %u\n", n_candidate);
#endif

	if (!n_candidate)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "No sequences to "
								"cluster.\n");

	if ((err = hdclust_malloc(opt, dat, ini, K_space, n_candidate)))
		return err;

	/* --------------------------------------------------------------------------- */
	/* Initialization: choose the most abundant unique sequences */

	/* select the first haplotype with the highest abundance */
	hd_update_seeds(dat, ini, K_curr, curr_useq_idx);

#ifdef DEBUG_AMPLICI
	debug_msg(DEBUG_I, fxn_debug, "Selecting %d (1/%d) with observed count"
		" %u\n", curr_useq_idx, n_candidate, ini->uniq_seq_count[i]);
	debug_call(DEBUG_II, fxn_debug, fprint_haplotype(stderr,
			ini->seeds[K_curr], ini->seed_lengths[K_curr]));
	debug_msg_cont(DEBUG_II, fxn_debug, "\n");
#endif

	ini->cluster_id[0] = K_curr++;

	/* ------------------------------------------------------------- */
	/* choose other haplotypes */

	while (curr_useq_idx < n_candidate) {

		unsigned int k_src = K_curr;

		/* If we need more space, we need to realloc the space */
		/* For increasing K */
		if (K_curr == K_space) {
#ifdef DEBUG_AMPLICI
			debug_msg(DEBUG_III, fxn_debug, "Begin reallocation");
#endif
			K_space = K_space + ini_K;
			if ((err = hdclust_realloc(opt, ini, mod, K_curr,
					K_space, dat->sample_size, n_candidates,
							dat->max_read_length)))
				return err;
#ifdef DEBUG_AMPLICI
			debug_msg(DEBUG_III, fxn_debug,
						"Finish reallocation\n");
#endif
		}

		/* select the haplotype temporarily */
		hd_update_seeds(dat, ini, K_curr, curr_useq_idx);

#ifdef DEBUG_AMPLICI
		debug_msg(DEBUG_I, fxn_debug, "Selecting %d/%d\n", curr_useq_idx,
								 n_candidate);
		debug_call(fxn_debug >= DEBUG_II, fxn_debug, fprint_haplotype(
			stderr, ini->seeds[K_curr], ini->seed_lengths[K_curr]));
		debug_msg_cont(DEBUG_II, fxn_debug, "\n");
#endif

		/* check if candidate is indel error of previous haplotype */
		if ((err = get_src_haplotype(opt, dat, mod, ini, K_curr + 1,
								&k_src)))
			return err;

		/* new haplotype not derived from previous haplotype */
		if (k_src == K_curr)
			++K_curr;
		ini->cluster_id[curr_useq_idx++] = k_src;
		++ini->cluster_size[k_src];

	}

	opt->K = K_curr;

	ini->cluster_membership = malloc(K_curr
					* sizeof(*ini->cluster_membership));
	if (!ini->cluster_membership)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
					"initializer::cluster_membership");

	for (unsigned int k = 0; k < K_curr; ++k) {
		ini->cluster_membership[k] = malloc(ini->cluster_size[k]
					* sizeof(**ini->cluster_membership));
		if (!ini->cluster_membership[k])
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"initializer::cluster_membership[%u]", k);
		ini->cluster_size[k] = 0;
	}

	for (unsigned int i = 0; i < n_candidates; ++i) {
		unsigned int k = ini->cluster_id[i];
		ini->cluster_membership[k][ini->cluster_size[k]++] = i;
	}

	opt->K_space = K_space;

	return err;
} /* cluster_by_hd */


/**
 * Allocate additional space specific for hd-clustering.
 *
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param ini		pointer to initializer object
 * @param array_fp	sequences of false positives
 * @param fp_abun	estiamted abundance of false positives
 * @param fp_trans	transition prob of false positive
 * @param nw_result	nw alignment result
 * @param nw_alen	length of nw alignemnt
 * @param K_space	space for K clusters
 * @param n_candidate	total number of candidates
 *
 * @return err		error status
 **/
int hdclust_malloc(options *opt, data *dat, initializer *ini,
				unsigned int K_space, unsigned n_candidate)
{
	/* haplotypes idx in unique sequence table
	 * need reallocation when K increases
	 */
	if (!ini->H)
		ini->H = malloc(K_space * sizeof(*ini->H));
	if (!ini->H)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "initializer::H");
	
	/* overwrite make_initializer() decision */
	if (ini->cluster_size) {
		mmessage(WARNING_MSG, NO_ERROR, "Not expecting to be freeing "
						"initializer::cluster_size.\n");
		free(ini->cluster_size);
	}

	ini->cluster_size = calloc(K_space, sizeof(*ini->cluster_size));
	if (!ini->cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::cluster_size");

	if (ini->cluster_id) {
		mmessage(WARNING_MSG, NO_ERROR, "Not expecting to be freeing "
						"initializer::cluster_id.\n");
		free(ini->cluster_id)
	}
	
	ini->cluster_id = malloc(K_space * sizeof(*ini->cluster_id));
	if (!ini->cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::cluster_id");

	if (ini->cluster_membership)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "initializer::"
			"cluster_membership should not be allocated.\n");
	
	return NO_ERROR;
} /* hdclust_malloc */


/**
 * Reallocate space when num of clusters increases above buffer.
 *
 * @param opt			pointer to options object
 * @param mod			pointer to model object
 * @param ini			pointer to initializer object
 * @param preK			current space for preK clusters
 * @param K			realloc space for K clusters
 * @param sample_size		number of total reads
 * @param hash_length		total number of candidates
 * @param max_read_length	length of reads
 *
 * @return			error status
 **/
int hdclust_realloc(options *opt, initializer *ini, model *mod,
	unsigned int preK, unsigned int K, size_t sample_size,
	unsigned int hash_length, unsigned int max_read_length)
{
	int err = NO_ERROR;

	unsigned int K_change = K - preK;

	if (K_change > 0) {

		/* H */
		unsigned int *H = realloc(ini->H, K * sizeof(*ini->H));

		if (!H)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"initializer::H");
		ini->H = H;

		unsigned int *cluster_size = realloc(ini->cluster_size,
						K * sizeof(*ini->cluster_size));

		if (!cluster_size)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::cluster_size");
		ini->cluster_size = cluster_size;

		unsigned int *cluster_id = realloc(ini->cluster_id,
						K * sizeof(*ini->cluster_id));

		if (!cluster_id)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer::cluster_id");

		ini->cluster_id = cluster_id;

		if ((err = realloc_seeds(ini, max_read_length, preK, K)))
			return err;

	}

	return err;

}/* hdclust_realloc */


/**
 * Decide whether a candidate haplotype should be included in the haplotype set.
 * - If \par options::nw_align == ALIGNMENT_HAPLOTYPES, then align the candidate
 *   to existing haplotypes and
 *   Check Hamming distance
 *
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param K		number of selected haplotypes, including candidate
 * @param k_src		index of source haplotype
 *
 * @return		error status
 *
 **/
int get_src_haplotype(options *opt, data *dat, model *mod, initializer *ini,
					unsigned int K, unsigned int *k_src)
{

#ifdef DEBUG_AMPLICI
	int fxn_debug = opt->info;
#endif
	int err = NO_ERROR;
	unsigned int curr_K = K - 1;

	*k_src = K;

	/* check for indel errors */
	if (opt->nw_align == ALIGNMENT_HAPLOTYPES) {
		if ((err = get_src_haplotype_with_indels(opt, dat, mod, ini,
								curr_K, k_src)))
			return err;
	} else if ((err = get_src_haplotype_without_indels(opt, dat, mod, ini,
								curr_K, k_src))) {
			return err;
	}
#if def D EBUG_AMPLICI
		debug_msg(DEBUG_I, fxn_debug, "indel error? %s.\n",
					k_src < K ? "yes" : "no");
#endif

	return err;
} /* get_src_haplotype */

/**
 * update seeds table in initializer when select a new haplotype
 *
 * @param ini		pointer to initializer object
 * @param dat		pointer to data object
 * @param K_prev	num of preselected haplotypes
 * @param idx		sequence index of the new haplotype
 *
 * @return err		error status
 *
 **/
int hd_update_seeds(data *dat, initializer *ini, unsigned int K_prev,
							 unsigned int idx)
{

	//ini->seed_idx[K_prev] = ini->uidx_to_ridx[idx]; // idx in dmat and qmat
	ini->seed_lengths[K_prev] = dat->lengths[ini->uidx_to_ridx[idx]];
	memcpy(ini->seeds[K_prev], dat->dmat[ini->uidx_to_ridx[idx]],
		dat->max_read_length * sizeof(**ini->seeds));
	ini->H[K_prev] = idx ; // idx in unique sequence table

	return NO_ERROR;
}/* hd_update_seeds */


/**
 * Check candidate haplotype for possible indel misread of existing haplotypes.
 * Specifically, align the candidate haplotype to all previous haplotypes and
 * compute edit distance.
 *
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param K_curr	number of haplotypes, excluding candidate
 * @param hap		source haplotype (maybe itself)
 *
 * @return err		error status
 **/
int get_src_haplotype_with_indels(options *opt, data *dat, model *mod, initializer *ini,
						unsigned int K_curr, int *hap)
{
	int err = NO_ERROR;
#ifdef DEBUG_AMPLICI
	int fxn_debug = opt->info;//ABSOLUTE_SILENCE;//DEBUG_IV;//DEBUG_III;//

	debug_msg(DEBUG_III, fxn_debug, "[scores] (ma=%i mm=%i gap=%i) band=%i "
		"sg=%i k=%i dbg=%i\n", opt->score[0][0], opt->score[0][1],
		opt->gap_p, opt->band, opt->ends_free, K_curr, fxn_debug);
#endif

	/* nw alignment for the candidate sequence and existing haplotypes */
	*hap = K_curr + 1;

	/* haplotype candidate */
	unsigned int rlen = ini->seed_lengths[K_curr];
	data_t *rseq = ini->seeds[K_curr];
	double min_ed = INFINITY;
	unsigned int min_k = K_curr;
	double scaling_const = 30.53628;		/* assume scores -3, -2, 2, -5 */
	double pgap = exp(opt->gap_p)/scaling_const;	/* assume iid scoring model */
	pgap = opt->indel_error;			/* change mind: overwrite with amplici assumption */

	/* align candidate haplotype to each existing haplotype */
	for (unsigned int k = 0; k < K_curr; k++) {
		double ascore = 0;
		size_t alen;
		unsigned int hap_len = ini->seed_lengths[k];
		unsigned int n_indels = 0, n_5prime = 0, n_gaps = 0, m_mismatch = 0;
		data_t *hap_seq = ini->seeds[k];

		data_t **aln = nwalign(hap_seq, rseq, (size_t) hap_len,
			(size_t) rlen, opt->score, opt->gap_p, opt->band,
				opt->ends_free, NULL, &err, &alen, &ascore);

#ifdef DEBUG_AMPLICI
		debug_call(fxn_debug >= DEBUG_V, fxn_debug,
					print_alignment(stderr, aln, alen));
#endif

		/* calculate number of indels and mismatch based on alignment */
		ana_alignment(aln, alen, rlen, &n_indels, &n_gaps, &n_5prime,
			&n_mismatch, opt->ends_free, ABSOLUTE_SILENCE);

#ifdef DEBUG_AMPLICI
		debug_msg(DEBUG_V, fxn_debug, "Haplotype %u score=%f alen=%zu "
				"id=%u idt=%u id5=%u mm=%u\n", k, ascore, alen,
					n_indels, n_gaps, n_5prime, n_mismatch);
#endif

		/* gap event = >0 insertions after haplotype position except
		 *     end or >0 deletions of consecutive haplotype nucleotides
		 * gap = insertion or deletion in haplotype or haplotype
		 * nw_indels[k] counts gap events
		 * n_indels counts gaps
		 * n_5prime counts gaps at 5' end
		 * be suspicious of alignment with more than 50% gap events,
		 * more than 2 clustered gaps unless they are at 5' end and we
		 * are ignoring these
		 */

		if ((n_indels + n_mismatch) / alen > min_ed) {
			min_ed = (n_indels + m_mismatch) / alen;
			min_k = k;
		} else if ((n_indels + n_mismatch) / alen < opt->hamming_distance) {
			debug_msg(1, 0, "h=%u also candidate\n", k);
		}
		if (aln) {
			if (aln[0])
				free(aln[0]);
			if (aln[1])
				free(aln[1]);
			free(aln);
		}
	}

	/*------------------------------------------------------------------- */
	/* recalculate transition prob for each read with same sequence as
	 * candidate haplotype now using NW alignment; need to allocate memory
	 */

	/* just a simple check for indel errors */
	if (min_ed <= opt->hamming_proportion)
		*hap = min_k;

	return err;
} /* get_src_haplotype_with_indels */


/**
 * Check candidate haplotype for possible misread of existing haplotypes.
 * Compute Hamming distance to all previous haplotypes and assign to closest
 * within options::hamming_proportion.
 *
 * @param opt		pointer to options object
 * @param dat		pointer to data object
 * @param mod		pointer to model object
 * @param ini		pointer to initializer object
 * @param K_curr	number of haplotypes, excluding candidate
 * @param hap		source haplotype (maybe itself)
 *
 * @return err		error status
 **/
int get_src_haplotype_with_indels(options *opt, data *dat, model *mod, initializer *ini,
						unsigned int K_curr, int *hap)
{
	int err = NO_ERROR;
#ifdef DEBUG_AMPLICI
	int fxn_debug = opt->info;//ABSOLUTE_SILENCE;//DEBUG_IV;//DEBUG_III;//

	debug_msg(DEBUG_III, fxn_debug, "hp=%f k=%i dbg=%i\n",
				opt->hamming_proportion, K_curr, fxn_debug);
#endif

	/* nw alignment for the candidate sequence and existing haplotypes */
	*hap = K_curr + 1;

	/* haplotype candidate */
	unsigned int rlen = ini->seed_lengths[K_curr];
	data_t *rseq = ini->seeds[K_curr];
	unsigned int min_k = 0;
	double min_hp = INFINITY

	/* align candidate haplotype to each existing haplotype */
	for (unsigned int k = 0; k < K_curr; k++) {
		unsigned int alen = MIN(ini->seed_length[k], rlen);
		data_t *hap_seq = ini->seeds[k];

		unsigned int hd = hamming_uchar_dis(hap_seq, rseq, alen);

#ifdef DEBUG_AMPLICI
		debug_msg(DEBUG_V, fxn_debug, "h=%u hd=%u hp=%f\n", k, hd,
								hd / alen);
#endif

		if (hd / alen > min_hp) {
			min_hp = hd / alen;
			min_k = k;
		} else if (hd / alen < opt->hamming_distance) {
			debug_msg(1, 0, "h=%u also candidate\n", k);
		}
		
	}

	if (min_hp <= opt->hamming_proportion)
		*hap = min_k;

	return err;
} /* get_src_haplotype_with_indels */


/* print reads assignment */
void fprint_assignment(FILE *fp, unsigned int *v, size_t n, unsigned int max, int width, int newline){
	size_t i;
	for (i = 0; i < n; ++i) {
		if (v[i]< max) {
			if (width)
				fprintf(fp, " %*u", width, v[i]);
			else
				fprintf(fp, " %u", v[i]);
		} else {
			fprintf(fp, " NA");
		}

	}
	if (newline) fprintf(fp, "\n");
} /* fprint_assignment */

void fprint_haplotype(FILE *fp, data_t *data, unsigned int len)
{
	for (size_t j = 0; j < len; ++j)
		fprintf(fp, "%c", xy_to_char[(int)data[j]]);
} /* fprint_haplotype */

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
			ri->optimal_cluster_id[i] = K; // treat those outliers as a new cluster

		ri->optimal_cluster_size[ri->optimal_cluster_id[i]]++;
	}

	return NO_ERROR;
}/* likelihood_filter */

/**
 * Reads assignment with given haplotype set
 *
 * @param ini	pointer to initializer object
 * @param dat	pointer data object
 * @param opt	pointer to options object
 * @param mod	pointer to model object
 * @param ri	pointer to run_info object
 *
 * return	err status
 **/
int reads_assignment(options * opt, data * dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;
#ifdef DEBUG_AMPLICI
	int fxn_debug = opt->info;
#endif
	//double l1third = 1./3;

	/* maybe use error profile */
	double *error_profile = NULL;
	if (opt->use_error_profile && mod->error_profile) {
		error_profile = mod->error_profile;
#ifdef DEBUG_AMPLICI
		debug_msg(DEBUG_II, fxn_debug, "Using error profile from "
					"'%s'. \n", opt->error_profile_name);
#endif
	}

	if ((err = trans_expectation(opt, mod, dat, ini, error_profile,
								mod->eik, 0)))
		return err;
	/* Keep codes below for debug purpose */
	/*

	for(unsigned int u = 0; u <dat->hash_length; ++u ){
		data_t *read = dat->dmat[ini->uidx_to_ridx[u]];
		unsigned int rlen = dat->lengths[ini->uidx_to_ridx[u]];

		unsigned int count = ini->uniq_seq_count[u]; // num. of reads wih unique seq
		size_t *idx_array; // idx of reads

		if ((err = find_index(dat->seq_count, read, rlen, &idx_array)))
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"Cannot find in the hash table !");

		// align to haplotypes
		for (unsigned int h = 0; h < opt->K; ++h) {

			data_t *hap_seq = ini->seeds[h];

			if (opt->nw_align == NO_ALIGNMENT) {

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
							double ep = dat->error_prob[dat->qmat[id][j]];
							if (dat->dmat[id][j] == hap_seq[j])
								eik += log(1 - ep);
							else
								eik += log(ep) + l1third;
						}
					}
					mod->eik[h*dat->sample_size+ id] = eik;
				}

			} else {
				size_t alen = dat->max_read_length;
				unsigned int nindels = 0;
				unsigned int nmismatch = 0;

				data_t **aln = nwalign(hap_seq, read,
				(size_t) ini->seed_lengths[h],
				(size_t) rlen,
				opt->score, opt->gap_p, opt->band, 1, NULL,
								&err, &alen, NULL);

				// count for number of indels
				ana_alignment(aln, alen, rlen, &nindels, NULL, NULL,
						&nmismatch, opt->info);

				for (unsigned int r = 0; r<count;++r) {

					mod->eik[h*dat->sample_size + idx_array[r]]
						= trans_nw(opt, mod, aln, alen,
						nmismatch, nindels,
						error_profile,
						mod->err_encoding,
						dat->qmat[idx_array[r]],
						dat->n_quality, rlen,
						dat->error_prob, opt->ends_free);

#ifdef DEBUG_AMPLICI
					debug_msg(DEBUG_III, fxn_debug, "num of indels: %i; num of "
							"mismatch: %i\n", nindels, nmismatch);
#endif

				}
				// free
				if (aln) {
					free(aln[0]);
					free(aln[1]);
					free(aln);
					aln = NULL;
				}
			}
		}
	}
	*/

	if (opt->trans_matrix) {
		FILE *fp = fopen(opt->trans_matrix, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->trans_matrix);
		//fprint_vectorized_matrix(fp, mod->eik, dat->sample_size, opt->K,0); not work
		for (size_t i = 0; i < dat->sample_size; ++i) {
			fprintf(fp, "%3lu", i);
			for (unsigned int j = 0; j < opt->K; ++j)
				fprintf(fp, " %8.2e", mod->eik[j*dat->sample_size + i]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	/* simply update mod->pi */
	assign_clusters(mod->eik, opt->K, dat->sample_size,
			ri->optimal_cluster_size, ri->optimal_cluster_id, 1);
	for (unsigned int k = 0; k < opt->K; ++k) {
		mod->pi[k] = (double) ri->optimal_cluster_size[k]
							/ dat->sample_size;
		if (!mod->pi[k])
			mod->pi[k] = 1.0 / dat->sample_size; // possible if given haplotypes
		mod->pi[k] = log(mod->pi[k]);
	}

	/* update mod->eik with new estimated mod->pi */
	for (unsigned int r = 0; r<dat->sample_size; ++r)
		for (unsigned int k = 0; k < opt->K; ++k)
			mod->eik[k*dat->sample_size+r] += mod->pi[k];

	/* reassign reads with updated mod->eik (unnormalized ) */
	assign_clusters(mod->eik, opt->K, dat->sample_size,
			ri->optimal_cluster_size, ri->optimal_cluster_id, 1);

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
