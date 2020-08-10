/**
 * @file ampliclust.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 * Cluster amplicon reads.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#include "ampliclust.h"
#include "initialize.h"
#include "statistics.h"
#include "lmath.h"
#include "io.h"
#include "hash.h"
#include "align.h"

#include "error.h"


/**
 * Assign reads to clusters based on posterior probabilities.
 *
 * @param eik	matrix E (contains posterior probabilities)
 * @param K 	number of clusters
 * @param n	number of reads
 * @param cs	cluster size (\par K x 1, to be written)
 * @param ci	cluster index (\par n x 1, to be writen)
 * @param by_K	matrix E is K x n (0) or n x K (1)
 * @return	error status
 */
int assign_clusters(double *eik, unsigned int K, size_t n,
		unsigned int *cs, unsigned int *ci, int by_K)
{
	double max;
	size_t i, k, l = 0;

	for (k = 0; k < K; ++k)
		cs[k] = 0;

	for (i = 0; i < n; ++i) {
		max = -INFINITY;
		for (k = 0; k < K; ++k)
			if ((by_K ? eik[n*k + i] : eik[K*i + k]) > max) {
				max = by_K ? eik[n*k + i] : eik[K*i + k];
				l = k;
			}
		/* three places to store result */
		ci[i] = l;
		++cs[l];
	}

	return NO_ERROR;
} /* assign_clusters */


/**
 * [KSD: stays here because uses data_t data format, which fastq.[ch] do not know]
 */
void fprint_fasta(FILE *fp, data_t *data, size_t n, size_t p, unsigned int*len, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		fprintf(fp, ">%s%lu\n", prefix, i);
		for (size_t j = 0; j < len[i]; ++j)
			fprintf(fp, "%c", xy_to_char[(int)data[i*p + j]]);
		fprintf(fp, "\n");
	}
} /* fprint_fasta */

/**
 * [KSD: see fprint_fasta()]
 */
#ifdef USE_CURSES
void wprint_fasta(WINDOW *wp, data_t *data, size_t n, size_t p, char const * const prefix) {
	for (size_t i = 0; i < n; ++i) {
		wprintw(wp, ">%s%u\n", prefix, i);
		for (size_t j = 0; j < p; ++j)
			wprintw(wp, "%c", xy_to_char[(int)data[i*p + j]]);
		wprintw(wp, "\n");
	}
} /* wprint_fasta */
#endif

/**
 * [KSD: stays here because uses data_t data format, which fastq.[ch] do not know]
 */
void fprint_alignment2(FILE *fp, data_t **data, size_t n, size_t p) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < p; ++j)
			fprintf(fp, "%c", xy_to_char[(int) data[i][j]]);
		fprintf(fp, "\n");
	}
} /* fprint_alignment2 */

/**
 * [KSD: see fprint_alignment2()]
 */
#ifdef USE_CURSES
void wprint_alignment2(WINDOW *wp, data_t **data, size_t n, size_t p) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < p; ++j)
			wprintw(wp, "%c", xy_to_char[(int) data[i][j]]);
		wprintw(wp, "\n");
	}
} /* wprint_alignment2 */
#endif

/**
 * Create run_info object.
 *
 * @param ri	pointer to run_info object
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 */
int make_run_info(run_info **ri, data *dat, options *opt)
{
	run_info *rio;
	*ri = malloc(sizeof **ri);
	unsigned int K = opt->K + 1;

	if (*ri == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "run_info");

	rio = *ri;

	/* initialize status variables*/
	rio->optimal_cluster_id = calloc(dat->sample_size,
		sizeof *rio->optimal_cluster_id);

	if (!rio->optimal_cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_cluster_id");

	rio->optimal_cluster_ll = calloc(dat->sample_size,
		sizeof *rio->optimal_cluster_ll);

	if (!rio->optimal_cluster_ll)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_cluster_ll");

	rio->optimal_cluster_size = calloc(K,
		sizeof *rio->optimal_cluster_size);

	if (!rio->optimal_cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"run_info::optimal_cluster_size");
	
	return NO_ERROR;
}/* make_run_info*/

/**
 * Reallocate run_info object for altered K or sample_size.
 *
 * @param ri		pointer to run_info object
 * @param sample_size	size of new sample
 * @param K		new K
 * @return		error status
 */
int realloc_run_info(run_info *ri, size_t sample_size, unsigned int K,
					int change_sample, int change_K)
{

	UNUSED(change_sample);
	UNUSED(sample_size);
	
	if (change_K) {

		unsigned int *optimal_cluster_size = realloc(
					ri->optimal_cluster_size,
				K * sizeof *ri->optimal_cluster_size);

		if (!optimal_cluster_size)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"realloc.run_info.optimal_cluster_size");

		ri->optimal_cluster_size = optimal_cluster_size;
	}

	return NO_ERROR;
}/* realloc_run_info */


/**
 * Delete run_info object
 */
void free_run_info(run_info *ri)
{
	if (ri){
		if(ri->optimal_cluster_id) free(ri->optimal_cluster_id);
		if(ri->optimal_cluster_size) free(ri->optimal_cluster_size);
		if(ri->optimal_cluster_ll) free(ri->optimal_cluster_ll);
		free(ri);
	}
}/*free_run_info*/
