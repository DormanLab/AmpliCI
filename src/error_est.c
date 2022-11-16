/**
 * @file error_est.c
 * @author Xiyu Peng
 *
 * Naive error estimation algorithm for AmpliCI
 *
 * TODO
 * - put the thresholds for detecting UMI collision under user control
 */

#include <limits.h>
#include "amplici.h"
#include "ampliclust.h"
#include "statistics.h"
#include "lmath.h"
#include "io.h"
#include "error.h"
#include "error_est.h"
#include "util.h"
#include "partition.h"
//#include <R.h>

#define DEBUG 0
#define WEIGHTS_LOESS
#define FILTER

int err_per_nuc(unsigned int nuc, unsigned int n_quality,
				unsigned int *count_sum, unsigned int *err_cnt,unsigned int self_lines[4],
				double *error_prob_temp, int if_loess, double *qual, double *error_profile);
int error_count_generator(options *opt,data *dat, model *mod, initializer *ini, run_info *ri);
int err_cnt_gen_wpartition(options *opt, data *dat, initializer *ini);

int error_profile_generator(options *opt, data *dat, model *mod,
					initializer *ini, run_info *ri)
{
	int err = NO_ERROR;

	int outputprofile = 1;
	int output_encoding = XY_ENCODING;

	/* number of distinct quality scores */
	unsigned int err_length = MAX_ASCII_QUALITY_SCORE
				- MIN_ASCII_QUALITY_SCORE + 1;
	unsigned int self_lines[4] = {0, 5, 10, 15}; // idx of lines for A->A, T->T, G->G, C->C

	double *error_profile  = NULL;  // For final output
	double *error_profile_p = NULL;  // estimated error rates

/* [KSD] What is the difference between model::n_quality and err_length?
 * [KSD] The former is the realized range in the given fastq file, the latter
 * [KSD] is the theoretical range given our hard-coded assumptions about FASTQ.
 * [KSD] Why use both?
 * [XY] Allow users to better interpret the output error profile ?
 */

	/* memory allocation for error profile */
	ini->err_cnt = calloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
				* mod->n_quality, sizeof *ini->err_cnt);

	error_profile = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
			* err_length * sizeof *error_profile);  // Used in the final output

	error_profile_p = malloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES  // only contains error rates estimated from loess fit
			* mod->n_quality * sizeof *error_profile_p);  // Used in the estimation step

	if (!error_profile || !error_profile_p || !ini->err_cnt)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"error profile generator ");

	/* find the optimum error count matrix for different K */
	/* place in ini->err_cnt */
	if (opt->partition_file) {
		if ((err = err_cnt_gen_wpartition(opt, dat, ini)))
			return err;
	} else if ((err = error_count_generator(opt, dat, mod, ini, ri)))
		return err;

	/* Below call regression function to predict errors */
	err = error_predict(ini->err_cnt, error_profile_p, dat->n_quality,
								self_lines);

	/* copy error_profile_p to error_profile and output */
	double one_third = 1.0/3.0;

	for (unsigned int nu1 = 0; nu1 < NUM_NUCLEOTIDES; nu1++) {
		for (unsigned int nu2 = 0; nu2 < NUM_NUCLEOTIDES; nu2++) {
			unsigned int r = nu1 * NUM_NUCLEOTIDES + nu2;

			for (size_t q = MIN_ASCII_QUALITY_SCORE;
				q <= MAX_ASCII_QUALITY_SCORE; ++q) {
				unsigned int idx = r * err_length
					+ q - MIN_ASCII_QUALITY_SCORE;

				if (q >= dat->fdata->min_quality
					&& q <= dat->fdata->max_quality) {
					error_profile[idx] = 1000 * exp(
						error_profile_p[
						r * dat->n_quality + q
						- dat->fdata->min_quality]);
				/* Use the literal error rate when there is no data */
				} else if (r == self_lines[nu1]) {
					error_profile[idx] = 1000 * (1 -
							error_prob(dat->fdata,
						(q - dat->fdata->min_quality)));
				} else {
					error_profile[idx] = 1000 * one_third
						* error_prob(dat->fdata,
						(q - dat->fdata->min_quality));
				}
			}
		}
	}

	if (outputprofile) {
		if (output_encoding == XY_ENCODING) {
			FILE *fp = fopen(opt->outfile_base, "w");

        		if (!fp)
            			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->outfile_base);

			fprint_error_profile(fp, error_profile, NUM_NUCLEOTIDES * NUM_NUCLEOTIDES, err_length);

			mmessage(INFO_MSG, NO_ERROR, "Output the error profile: %s \n", opt->outfile_base);

	       		fclose(fp);
		} else {
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "Codes are not ready");
		}
	}

	free(error_profile);
	free(error_profile_p);

	return err;
}/* error_profile_generator */

/**
 * iteratively generate error count profile, which is stored in ini->err_cnt
 *
 * @param opt 	options object
 * @param dat	data object
 * @param mod	error model object
 * @param ini	initializer object
 * @param ri	run initializer
 */
int error_count_generator(options *opt, data *dat, model *mod,
					initializer *ini, run_info *ri)
{

	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;
	int fix_K = opt->K_fix_err;
	unsigned int K_max;
	unsigned int nrow = NUM_NUCLEOTIDES * NUM_NUCLEOTIDES;
	unsigned int K_space = opt->K;

	unsigned int self_lines[4] = {0, 5, 10, 15};

	/* maximum K used for error estimation */
	if (fix_K)
		K_max = opt->K;	// start with up to options::K haplotypes
	else
		K_max = 1;	// start with K = 1 haplotypes

	/* allocate memory for err_cnt */
	unsigned int *err_cnt_prev = NULL;
	unsigned int *err_cnt_cur = NULL;

	err_cnt_prev = calloc(nrow * mod->n_quality, sizeof *err_cnt_prev);

	if (!err_cnt_prev)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "error_cnt_prev");

	err_cnt_cur = calloc(nrow * mod->n_quality, sizeof *err_cnt_cur);

	if (!err_cnt_cur) {
		free(err_cnt_prev);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "error_cnt_cur");
	}

	/* for error profile */
	if (!mod->error_profile)
		mod->error_profile = calloc(nrow * mod->n_quality,
						sizeof *mod->error_profile);
	if (!mod->error_profile) {
		free(err_cnt_prev);
		free(err_cnt_cur);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "error_profile");
	}

	do {

		/* initialize parameters for rerunning AmpliCI */
		opt->K = K_space;

		for (unsigned int i = 0; i < nrow; i ++)
			for (unsigned int q = 0; q < mod->n_quality; ++q)
				err_cnt_cur[i*mod->n_quality + q] = 0;

		/* select up to K_max haplotypes */
		if ((err = haplotype_selection(opt, dat, mod, ini, K_max)))
			return err;

		if ((err = realloc_run_info(ri, dat->sample_size, opt->K + 1,
									0, 1)))
			return err;

		/* assign reads to haplotypes using posterior assignment
		 * probability
		 */
		assign_clusters(mod->eik, opt->K, dat->sample_size,
			ri->optimal_cluster_size, ri->optimal_cluster_id, 1);

		/* count number of reads that have maximum posterior probability of assignment > cutoff */
		unsigned int count = 0;
		for (unsigned int i = 0 ; i < dat->sample_size; i++){
			ri->optimal_cluster_ll[i] =
				mod->pi[ri->optimal_cluster_id[i]]
					+ ini->e_trans[ri->optimal_cluster_id[i]
							* dat->sample_size + i];
			if (ri->optimal_cluster_ll[i] > opt->ll_cutoff)
				++count;
		}

		mmessage(INFO_MSG, NO_ERROR, "number of reads that have "
						"ll> ll_cutoff:%i \n", count);

		/* use estimated errors for next haplotype addition */
		/* [TODO] It is arbitrary to set the cutoff is 0.5 */

		if (count > dat->sample_size * 0.5)
			opt->use_error_profile = 1;

		/* likelihood filter: filtered into group options::K */
		#ifdef FILTER
		if(opt->filter_reads)
			likelihood_filter(opt->K, opt->ll_cutoff, NULL, mod->pi, ini->e_trans,
							dat->sample_size, ri);
		#endif

		/*  Generate current error count profile.
		 *	count all (n,m,q) combinations for true base n, read base m
		 * and quality score q, assuming assigned haplotype is true
		 */
		err_count_with_assignment(dat, ri->optimal_cluster_id, ini->seeds,ini->seed_lengths, err_cnt_cur,opt->K);
		/*
		for (size_t i = 0; i < dat->sample_size; ++i) {
			// ignore any filtered
			if (ri->optimal_cluster_id[i] == opt->K)
				continue;

			for (size_t j = 0; j < dat->lengths[i]; ++j)
				err_cnt_cur[dat->n_quality * (
					ini->seeds[ri->optimal_cluster_id[i]][j]
					* NUM_NUCLEOTIDES + dat->dmat[i][j])
					+ dat->qmat[i][j]]++;// order of A,C,T,G
		} */

		if (K_max == opt->K_max)
			break;

		/* compare the distance between err_cnt_cur and err_cnt_pre */
		double min_cos_dist = cos_dist(err_cnt_prev, err_cnt_cur, nrow,
								mod->n_quality);

		/* cosine distance close enough to 1 terminates the loop */
		if (min_cos_dist > opt->min_cosdist) {

			mmessage(INFO_MSG, NO_ERROR, "Convergence at K = %i "
				"with minimum distance %3.8f\n", opt->K,
							exp(min_cos_dist));

			break;
		}

		/* new error profile */
		if ((err = error_predict(err_cnt_cur, mod->error_profile,
						mod->n_quality, self_lines)))
			return err;

		/* copy err_cnt_cut to err_cnt_prev */
		memcpy(err_cnt_prev, err_cnt_cur,
				nrow * mod->n_quality * sizeof(*err_cnt_prev));
		++K_max;

	} while (!fix_K);

	/* output the result we are interested */
	if (fxn_debug > DEBUG_I) {
		FILE *fp = fopen(opt->outfile_base, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->outfile_base);

		fprintf(fp, "error count:");
		fprint_vectorized_uintmatrix(fp, err_cnt_cur,
			NUM_NUCLEOTIDES * NUM_NUCLEOTIDES, mod->n_quality, 1);

		fclose(fp);
	}

	/* copy the final err_count to ini->err_count */
	memcpy(ini->err_cnt, err_cnt_cur, nrow * mod->n_quality
						* sizeof(*ini->err_cnt));

	if (err_cnt_prev)
		free(err_cnt_prev);
	if (err_cnt_cur)
		free(err_cnt_cur);

	return err;
}/* error_count_generator */

/**
 * Calculate the minimum cosine distance (log version) between two sets of
 * error counts.
 *
 * @param mat1	count of previous error and quality score combinations
 * @param mat2	count of new error and quality score combinations
 * @param nrow	number of rows in matrices
 * @param ncol	number of columns in matrices
 */
double cos_dist(unsigned int *mat1, unsigned int *mat2, unsigned int nrow,
								unsigned int ncol)
{

	double epsilon = 1e-50;
	double min_cos_dist = 0.;
	unsigned int sum1 = 0, sum2 = 0, idx;

	double *cos_dist = calloc(ncol, sizeof *cos_dist);
	unsigned int *colsum1 = calloc(ncol, sizeof *colsum1);
	unsigned int *colsum2 = calloc(ncol, sizeof *colsum2);

	if (!cos_dist || !colsum1 || !colsum2) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "cosine distance");
		min_cos_dist = log(epsilon);	/* [KSD] Communicate lethal problem! */
		goto EXIT_COS_DIST;
	}

	/* column sum of the two matrix (avoid possible overflow) */
	/* [KSD] How do you avoid overflow? What are sum1 and sum2 doing? */
	for (unsigned int j = 0; j < ncol; j++) {
		colsum1[j] = 0;
		colsum2[j] = 0;
		for (unsigned int i = 0; i < nrow; i ++) {
			idx = i * ncol + j;
			colsum1[j] += mat1[idx];
			colsum2[j] += mat2[idx];
			if (!sum1 || !sum2) {
				sum1 += mat1[idx];
				sum2 += mat2[idx];
			}
		}
	}

	if (!sum1 || !sum2){
		min_cos_dist = log(epsilon);
		goto EXIT_COS_DIST;
	}

	/* cos distance */
	for (unsigned int j = 0; j < ncol; j++) {
		cos_dist[j] = 0.;
		double nu = 0.0, de1 = 0.0,de2 = 0.0;
		for (unsigned int i = 0; i < nrow; i++) {
			idx = i * ncol + j;
			/* avoid overflow of unsigned int multiplication */
			nu += ((double) mat1[idx] / colsum1[j])
					* ((double) mat2[idx] / colsum2[j]);
			de1 += ((double) mat1[idx] / colsum1[j])
					* ((double) mat1[idx] / colsum1[j]);
			de2 += ((double) mat2[idx] / colsum2[j])
					* ((double) mat2[idx] / colsum2[j]);
		}
		cos_dist[j] = log(nu + epsilon)
					- log(sqrt(de1 * de2) + epsilon);
		/* use the minimum value */
		if (cos_dist[j] < min_cos_dist)
			min_cos_dist = cos_dist[j];
	}

EXIT_COS_DIST:
	if (cos_dist)
		free(cos_dist);
	if (colsum1)
		free(colsum1);
	if (colsum2)
		free(colsum2);

	return min_cos_dist;
}/* cos_dist */

/* print error profile   */
void fprint_error_profile(FILE *fp, double *mat, unsigned int n, unsigned int l){
	unsigned int i, j;

	for (i = 0; i < n; ++i)
		for (j = 0; j < l; ++j)
			fprintf(fp, " %8.8f,", mat[i*l + j]);
} /* fprint_error_profile */

/**
 * Predict error rates based on the observed error counts.
 *
 * @param err_cnt	     error counts
 * @param error_profile	 error profile to be computed
 * @param n_quality	     number of distinct quality scores
 * @param self_lines	 idx of lines of self-transitions, like A->A, T->T, C->C, G->G
 * @return		         error status
 */
int error_predict(unsigned int *err_cnt, double *error_profile,
			unsigned int n_quality, unsigned int self_lines[4])
{

	int err = NO_ERROR;

	double *error_prob_temp = NULL;
	unsigned int *count_sum = NULL;
	double *qual = NULL;

	int if_loess = 1;

	/* malloc work space */
	error_prob_temp = malloc(n_quality * sizeof *error_prob_temp);
	count_sum = calloc(n_quality, sizeof *count_sum);    // initial value: 0
	qual = calloc(n_quality, sizeof *qual);

	if (!error_prob_temp || !count_sum || !qual) {
		err = mmessage(ERROR_MSG, MEMORY_ALLOCATION, "error_predict");
		goto EXIT_ERROR_PREDICT;
	}

	for (unsigned int q = 0; q < n_quality; q++)
		qual[q] = (double) q;

	/* estimate error profile for A, C, T, G each nuc */
	for (unsigned int nuc = 0; nuc < NUM_NUCLEOTIDES; nuc ++)
		if ((err = err_per_nuc(nuc, n_quality, count_sum, err_cnt,
				self_lines, error_prob_temp, if_loess, qual,
								error_profile)))
			return err;

EXIT_ERROR_PREDICT:
	if (error_prob_temp)
		free(error_prob_temp);
	if (count_sum)
		free(count_sum);
	if (qual)
		free(qual);

	return NO_ERROR;
}/* error_predict */

/* function to call loess regression */
int loess_est(loess *lo,double *y,double *x, unsigned int* weights, unsigned int n){

	double span = 0.75;

	loess_setup(x,y,(int)n,1,lo);
	lo->model.span = span;

#ifdef WEIGHTS_LOESS
	/* Use weights */
	for(unsigned int i = 0; i < n; i ++)
		lo->inputs.weights[i] = (double) weights[i] + 0.001; //avoid the problems
#endif

	loess_fit(lo);

	// check for bugs
#ifndef DEBUG
	loess_summary(lo);
#endif

	return NO_ERROR;
}/* loess_est */

/**
 * Estimate error profile for nucleotides A, T, C, G
 *
 * @param nuc			current nucleotide
 * @param n_quality		number of distinct quality scores
 * @param count_sum		room for one column of err_cnt
 * @param err_cnt		counts of all (n,m,q) combinations
 * @param self_lines		idx of lines of self-transitions, like A->A, T->T, C->C, G->G
 * @param error_prob_temp	vector of length n_quality
 * @param if_loess		use loess regression
 * @param qual			all possible quality values
 * @param error_profile		estimated error profile (to compute)
 * @return			error status
 */
int err_per_nuc(unsigned int nuc, unsigned int n_quality,
			unsigned int *count_sum, unsigned int *err_cnt,
	unsigned int self_lines[4], double *error_prob_temp, int if_loess,
					double *qual, double *error_profile)
{
	UNUSED(if_loess);

	/* structs for regression */
	loess *lo = NULL;
	prediction *pre = NULL;

	unsigned int lb = NUM_NUCLEOTIDES * nuc;
	unsigned int ub = NUM_NUCLEOTIDES * nuc + NUM_NUCLEOTIDES;
	unsigned int idx;
	double lquarter = log(0.25);
	double lepsilon = log(1e-7);
	double l10 = log(10);

	/* total true nucleotides nuc observed with each quality score */
	for (unsigned int q = 0; q < n_quality; q++) {
		count_sum[q] = 0;
		for (unsigned int k = lb; k < ub; k++)
			count_sum[q] += err_cnt[k * n_quality + q];
	}

	for (unsigned int k = lb; k < ub; k++) {

		if (k == self_lines[nuc])      // skip self-transitions
			continue;

		/* prob. of each quality score when true nucleotide is nuc */
		for (unsigned int q = 0; q < n_quality; ++q)
			error_prob_temp[q] = log10((double)
				(err_cnt[k*n_quality + q] + 1)
				/ (count_sum[q] + 4));   //avoid numerical bugs

		/* malloc space for lo and pre */
		lo = malloc(sizeof *lo);
		if (!lo)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"loess object");

		pre = malloc(sizeof *pre);
		if (!pre) {
			free(lo);
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"prediction object");
		}

		/* weighted loess regression */
		loess_est(lo, error_prob_temp, qual, count_sum,
							n_quality);

		predict(qual, n_quality, lo, pre, 0);

		for (unsigned int q = 0; q < n_quality; ++q) {
			idx = k * n_quality + q;
			error_profile[idx] = pre->fit[q] * l10;  //log(error_prob)
			if (error_profile[idx] > lquarter)
				error_profile[idx] = lquarter; //avoid bugs
			if (error_profile[idx] < lepsilon)
				error_profile[idx] = lepsilon;
		}


		if (lo)
			loess_free_mem(lo);
		if (pre)
			pred_free_mem(pre);

		lo = NULL;
		pre = NULL;

	}

	/* prob of non-error for each quality score */
	for (unsigned int q = 0; q < n_quality; q++) {
		idx = self_lines[nuc] * n_quality + q;
		error_profile[idx] = 1.0;
		for (unsigned int k = lb; k < ub; k++)
			if (k != self_lines[nuc])
				error_profile[idx] -= exp(error_profile[k*n_quality+q]);
		error_profile[idx] = log(error_profile[idx]); // log version
	}

	return NO_ERROR;
}/* err_per_nuc */

/**
 * Generate error count profile (ini->err_cnt) with given partition.
 *
 * @param opt 	options object
 * @param dat	data object
 * @param mod	error model object
 * @param ini	initializer object
 * @param ri	run initializer
 * @return	error status
 */
int err_cnt_gen_wpartition(options *opt, data *dat, initializer *ini)
{
	int err = NO_ERROR;

	/*  read the partition file, start from 0 */
	if ((err = read_partition_file(opt->partition_file, ini->cluster_id,
							dat->sample_size)))
		return err;

	unsigned int max = 0;

	for (unsigned int i = 0; i < dat->sample_size; i++)
		if (ini->cluster_id[i] > max)
			max = ini->cluster_id[i];

	unsigned int K_partitions = max + 1;  // new K, number of UMI clusters

	/* find the most abundant sequence for each partition */
	hash **hash_list = NULL;
	hash_list = malloc(K_partitions * sizeof(*hash_list));
	if (!hash_list)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"Err_cnt.hash_list");

	/* use existing array to store number of haplotypes in each UMI subset */
	ini->cluster_size = malloc(K_partitions * sizeof(*ini->cluster_size));

	if (!ini->cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer:cluster_size");

	for (unsigned int k = 0; k < K_partitions; k++) {
		hash_list[k] = NULL;
		ini->cluster_size[k] = 0;
	}

	// check bugs in hash (currently not find any bugs )
	// unsigned int slen;
	for (size_t i = 0; i < dat->sample_size; ++i) {
		add_sequence(&hash_list[ini->cluster_id[i]],
				dat->dmat[i], dat->lengths[i], i, &err);
		if (err)
			return err;
	}

	// sort by abundance within each partition
	for (unsigned int k = 0; k < K_partitions; k++)
		sort_by_count(&hash_list[k]);

	// simple idea to detect collision in each UMI clusters
	// consider every member with abundance at least half of most
	// abundant and at least 2 to be a true molecule
	// [TODO] Put these choices under user control!
	// remember need to update both ini->cluster_id, ini->seeds

	// determine a new K

	unsigned int K_seeds = 0;

	for (unsigned int k = 0; k < K_partitions; ++k) {
		if (!hash_list[k])
			continue;

		unsigned int max_abun = hash_list[k]->count;  // the most abundant sequence
		unsigned int thres = (max_abun + 1) / 2;

		for (hash *s = hash_list[k]; s != NULL; s = s->hh.next) {
			if (s->count >= thres && s->count
					>= opt->seed_min_observed_abundance) {
				++K_seeds;
				++ini->cluster_size[k];

				// if we just need the most abundant sequence for each UMI cluster
				if (!opt->umicollision)
					break;

			} else {
				break;
			}
			//
			//if(s->count > max_abun){
			//	if(s->count > 20) fprintf(stderr, "s->count > 20\n");
			//	max_abun = s->count;
			//	idx = s->idx;
			//	memcpy(ini->seeds[k], s->sequence,
			//	dat->max_read_length * sizeof **ini->seeds);
			//}
		}
	}

	//unsigned int sum_subK;
	//for (unsigned int k =0 ; k < K_partitions; k++){
	//	sum_subK +=  ini->cluster_size[k];
	//}
	//fprintf(stderr, "sum_subK : %d\n", sum_subK);

	/* some clusters have no sufficiently replicate haplotype to be considered candidate */
	unsigned int num_null = 0;
	for (unsigned int k = 0; k < K_partitions; ++k)
		if (ini->cluster_size[k] == 0)
			num_null++;

	K_seeds += (opt->exclude_low_abundance_seeds ? 0 : num_null);

	if (!K_seeds)
		return mmessage(ERROR_MSG, NO_DATA, "No subset of UMI partition"
			" contains a haplotype observed at least %u times.  "
			"Cannot estimate error profile!\n",
			opt->seed_min_observed_abundance);

	/* realloc space for seeds */
	if (K_seeds != opt->K)
		if ((err = realloc_seeds(ini, dat->max_read_length,
							opt->K, K_seeds)))
			return err;

	/* identify true assignment and seeds when recognizing possible collision */

	unsigned int current_k = 0;
	unsigned int kept_reads = 0;
	unsigned int extra_k = 0;
	unsigned int *new_cluster_id = NULL;

	if (opt->exclude_low_abundance_seeds) {
		new_cluster_id = malloc(dat->sample_size
						* sizeof(*new_cluster_id));
		if (!new_cluster_id)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"new_cluster_id");
	}

	for (unsigned int k = 0; k < K_partitions; k++) {

		if (!hash_list[k])
			continue;

		hash *s = hash_list[k];
		unsigned int subK = ini->cluster_size[k];

		// no collision or zero cluster
		if (subK == 1 || (subK == 0
				&& !opt->exclude_low_abundance_seeds)) {

			memcpy(ini->seeds[current_k], s->sequence,
				dat->lengths[s->idx] * sizeof **ini->seeds);	
			ini->seed_lengths[current_k] = dat->lengths[s->idx];
			/* some clusters are dropped */
			if (opt->exclude_low_abundance_seeds)
				for (unsigned int i = 0; i < dat->sample_size;
									++i)
					if (ini->cluster_id[i] == k)
						new_cluster_id[i] = current_k;
			++current_k;

		/* exclude clusters without obvious centers */
		} else if (subK == 0 && opt->exclude_low_abundance_seeds) {
			for (unsigned int i = 0; i < dat->sample_size; ++i)
				if (ini->cluster_id[i] == k)
					new_cluster_id[i] = K_seeds;
		} else {  // collision
			unsigned int lidx[subK];

			// copy seeds information
			for (unsigned int sk = 0; sk < subK; ++sk) {
				//fprintf(stderr, "count: %d\n",s->count);

				if (sk == 0) {
					memcpy(ini->seeds[current_k],
								s->sequence,
							dat->lengths[s->idx] 	
							* sizeof(**ini->seeds));
					ini->seed_lengths[current_k]
							= dat->lengths[s->idx];
					lidx[sk] = current_k;
					++current_k;
				} else {
					unsigned int lk =
						opt->exclude_low_abundance_seeds
						? current_k++
						: K_partitions + extra_k++;
					//fprintf(stderr, "%d, %d\n",pos,K_seeds);
					memcpy(ini->seeds[lk], s->sequence,
							dat->lengths[s->idx]  	
							* sizeof(**ini->seeds));
					ini->seed_lengths[lk]
						= dat->lengths[s->idx];
					lidx[sk] = lk;
				}
				kept_reads += s->count;

				s = s->hh.next;
			//	fprintf(stderr, "%d, %d\n",pos,subK);
			}
			//fprintf(stderr, "%d, %d\n",add_num,subK-1);

			// assign members of split clusters by Hamming distance
			unsigned int min_dist;

			for (unsigned int i = 0; i < dat->sample_size; ++i) {

				if (ini->cluster_id[i] != k)
					continue;

				// [TODO] find the minmimal length of seeds and use it when calculate the hamming distance.
				min_dist = UINT_MAX;
				for (unsigned int sk = 0; sk < subK; sk++) {
					unsigned char *hapk
							= ini->seeds[lidx[sk]];
					unsigned int hdist = hamming_uchar_dis(
							dat->dmat[i], hapk,
							dat->max_read_length);	/* [BUG,KSD] Looks like a bug for variable length reads */

					// fprintf(stderr, "%d, %d\n",lidx[sk], hdist);
					if (hdist < min_dist) {
						if (opt->exclude_low_abundance_seeds)
							new_cluster_id[i] = lidx[sk];
						else
							ini->cluster_id[i] = lidx[sk];
						min_dist = hdist;
					}
				}
				//fprintf(stderr, "previous %d, current %d\n",k, ini->cluster_id[i]);
			}
		}
	}

	if (opt->exclude_low_abundance_seeds && current_k != K_seeds)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "current_k = %u != K_seeds = %u\n", current_k, K_seeds);
	if (opt->exclude_low_abundance_seeds) {
		free(ini->cluster_id);
		ini->cluster_id = new_cluster_id;
	}

	mmessage(INFO_MSG, NO_ERROR, "Will estimate errors from %u reads in "
		"%u out of %u clusters.\n", kept_reads, K_seeds, K_partitions);

	//for(unsigned int i = 0; i < dat->sample_size; i++){
	//	fprintf(stderr, "%d",ini->cluster_id[i]);
	//}


	/*
	for (unsigned int k = 0; k < K_partitions; k++){
		unsigned int max_abun = 0;
		hash *s = NULL;
		unsigned int idx = 0;
		// fprintf(stderr, "count : %d\n", hash_list[k]->count);
		for (s = hash_list[k]; s != NULL; s = s->hh.next) {
			//if(k==3)
			//	fprintf(stderr, "s->count : %d for the %d th \n", s->count, k);
			if(s->count > max_abun){
				//if(s->count > 20) fprintf(stderr, "s->count > 20\n");
				max_abun = s->count;
				idx = s->idx;
				//memcpy(ini->seeds[k], s->sequence,
				//	dat->max_read_length * sizeof **ini->seeds);
			}
		}
		memcpy(ini->seeds[k], dat->dmat[idx],
				dat->max_read_length * sizeof **ini->seeds);
		ini->seed_lengths[k] = dat->lengths[idx];
		if(hash_list[k])
			delete_all(&hash_list[k]);
	} */


	/* count error types with given assignment */
	err_count_with_assignment(dat, ini->cluster_id, ini->seeds,
				ini->seed_lengths, ini->err_cnt, K_seeds);

	/* free space */
	for (unsigned int k = 0; k < K_partitions; k++)
		if (hash_list[k])
			delete_all(&hash_list[k]);
	if (hash_list)
		free(hash_list);

	return NO_ERROR;
}/* err_cnt_gen_wpartition */


/**
 * count error types with given assignment
 *
 * @param dat            data object
 * @param cluster_id     reads assignment id
 * @param seeds          haplotypes
 * @param count_mat      array for count
 * @param K              maximum number of clusters
 * @return               error status
 */
int err_count_with_assignment(data *dat, unsigned int *cluster_id,
	data_t** seeds, unsigned int *seed_lengths, unsigned int *count_mat,
								unsigned int K)
{
	unsigned int kept_errors = 0, kept_bases = 0;

	for (size_t i = 0; i < dat->sample_size; ++i) {
		/* ignore any filtered reads */
		if (cluster_id[i] >= K)
			continue;

		for (size_t j = 0; j < dat->lengths[i]; ++j)
			if(j < seed_lengths[cluster_id[i]]) {
				count_mat[dat->n_quality * (
					seeds[cluster_id[i]][j]
					* NUM_NUCLEOTIDES + dat->dmat[i][j])
					+ dat->qmat[i][j]]++;// order of A,C,T,G
				if (seeds[cluster_id[i]][j] != dat->dmat[i][j])
					++kept_errors;
				++kept_bases;
		}
	}

	mmessage(INFO_MSG, NO_ERROR, "Estimating with %u observed errors out "
			"of %u kept base calls.\n", kept_errors, kept_bases);
	return NO_ERROR;
}/* err_count_with_assignment */

