/**
 * @file data.c
 * @author Karin S. Dorman
 *
 * Manipulate data structure for ampliclust.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "options.h"
#include "model.h"
#include "util.h"
#include "lmath.h"
#include "error.h"
#include "hash.h"

int build_hash(data *dat);

/**
 * Setup data object.
 */
int make_data(data **dat, options *opt)
{
	UNUSED(opt);
	data *dp;
	*dat = malloc(sizeof **dat);

	if (*dat == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data object");

	dp = *dat;

	dp->fdata = NULL;
	dp->coverage = NULL;
	dp->n_quality = 0;
	dp->dmat = NULL;
	dp->qmat = NULL;

	dp->error_prob = NULL;

	dp->read_idx = NULL;
	dp->sample_size = 0;
	dp->max_read_length = 0;
	dp->min_read_length = 0;
	dp->lengths = NULL;

	dp->offset = NULL;
	dp->max_offset = 0;

	dp->max_read_position = 0;

	//dp->true_cluster_id = NULL;
	//dp->true_cluster_size = NULL;
	
	dp->seq_count = NULL;
	dp->hash_length = 0;

	return NO_ERROR;
} /* make_data */

/**
 * Sync state of data object.
 */
int sync_state(data *dat, options *opt)
{
	int fxn_debug = ABSOLUTE_SILENCE;
	int err = NO_ERROR;
	UNUSED(opt);

	/* currently, the sample is the complete data */
	dat->sample_size = dat->fdata->n_reads;

	dat->max_read_length = dat->fdata->n_max_length;
	dat->min_read_length = dat->fdata->n_min_length;
	dat->max_read_position = dat->max_read_length;
	/* warning message */
	if (dat->max_read_length != dat->min_read_length)
		mmessage(WARNING_MSG, NO_ERROR,
			"Unequal reads length! \n");

	dat->n_quality = dat->fdata->max_quality - dat->fdata->min_quality + 1;

	dat->error_prob = malloc(dat->n_quality * sizeof * dat->error_prob);
	if(!dat->error_prob)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,"dat.error_prob");

	for (unsigned int q = 0; q < dat->n_quality; q++)
		dat->error_prob[q] = error_prob(dat->fdata, q);

	dat->lengths = malloc(dat->sample_size * sizeof *dat->lengths);

	if (!dat->lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.lengths");

	/* [BUG?] first sample are first data::sample_size reads, assumption
	 * continues throughout this function with consequences
	 */
	if (dat->fdata->n_lengths)
		memcpy(dat->lengths, dat->fdata->n_lengths, dat->sample_size
			* sizeof *dat->lengths);
	else
		for (size_t i = 0; i < dat->sample_size; ++i)
			dat->lengths[i] = dat->max_read_length;

	
	/* allocate the index array of reads */
	dat->read_idx = malloc(dat->fdata->n_reads * sizeof *dat->read_idx);

	if (!dat->read_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.read_idx");

	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		dat->read_idx[i] = i;

	/* allocate short-cut pointers to nucleotide sequences */
	dat->dmat = malloc(dat->sample_size * sizeof *dat->dmat);

	if (!dat->dmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.dmat");

	unsigned char *rptr = dat->fdata->reads;
	for (size_t i = 0; i < dat->sample_size; ++i) {
		dat->dmat[i] = rptr;
		rptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) sequence matrix\n",
		  dat->sample_size, dat->max_read_length);

	/* allocate short-cut pointers to quality sequences  */
	dat->qmat = malloc(dat->sample_size * sizeof *dat->qmat);

	if (!dat->qmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.qmat");

	unsigned char *qptr = dat->fdata->quals;
	for (size_t i = 0; i < dat->sample_size; i++) {
		dat->qmat[i] = qptr;
		qptr += dat->lengths[i];
	}

	debug_msg(DEBUG_I, fxn_debug, "Allocated %dx(%d) quality matrix\n",
		  dat->sample_size, dat->max_read_length);

	/*
	if (opt->offset_file) {
		dat->offset = malloc(dat->fdata->n_reads * sizeof *dat->offset);
		if (!dat->offset)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"data.offsets");
		if ((err = read_uints(opt->offset_file, dat->offset,
			dat->fdata->n_reads)))
			return err;
		for (size_t i = 0; i < dat->fdata->n_reads; ++i) {
			if (dat->offset[i] > dat->max_offset)
				dat->max_offset = dat->offset[i];
			if (dat->lengths[i] + dat->offset[i]
				> dat->max_read_position)
				dat->max_read_position = dat->lengths[i]
					+ dat->offset[i];
		}

	}

	debug_msg(DEBUG_I, fxn_debug, "dat->max_offset:%i\n", dat->max_offset);

	dat->coverage = malloc(dat->max_read_position * sizeof *dat->coverage);

	if (!dat->coverage)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.coverage");

	if(opt->offset_file){
		for (size_t j = 0; j < dat->max_read_position; ++j)
			dat->coverage[j] = 0;

		for (size_t i = 0; i < dat->sample_size; ++i)
			for (size_t j = 0; j < dat->lengths[i]; ++j)
				dat->coverage[j + (dat->offset ? dat->offset[i] : 0)]++;

		debug_msg(DEBUG_I, fxn_debug, "coverage");
	}else{
		for (size_t j = 0; j < dat->max_read_position; ++j)
			dat->coverage[j] = dat->sample_size;
	}
	*/
	
	if (!dat->fdata->empty)
		err = sync_data(dat, opt);

	return err;
} /* sync_state */


/**
 * Build hash of reads.
 *
 * @param dat	pointer to data object
 * @return	error status
 */
int build_hash(data *dat)
{
	int err = NO_ERROR;
	/* build hash table */
	dat->hash_length = 0;
	for (size_t i = 0; i < dat->sample_size; ++i) {
//fprintf(stderr, "%zu: %s\n", i, display_sequence(dat->dmat[i], dat->lengths[i], dat->fdata->read_encoding));
		dat->hash_length += add_sequence(&dat->seq_count, dat->dmat[i],
							dat->lengths[i], i);
	}

	/* store index of reads for all unique sequences */
	for (size_t i = 0; i < dat->sample_size; ++i)
		if ((err = add_seq_idx(dat->seq_count, dat->dmat[i],
							dat->lengths[i], i)))
			return err;
			
	/* sort hash table by count */
	sort_by_count(&dat->seq_count);

	return err;
} /* build_hash */


/**
 * Free data object.
 */
void free_data(data *dat)
{
	if (dat) {
		if (dat->fdata) free_fastq(dat->fdata);
		if (dat->dmat) free(dat->dmat);
		if (dat->qmat) free(dat->qmat);
		if (dat->read_idx) free(dat->read_idx);
		if (dat->lengths) free(dat->lengths);
		if (dat->offset) free (dat->offset);
		if (dat->coverage) free (dat->coverage);
		if (dat->seq_count) delete_all(dat->seq_count);
		if (dat->error_prob) free(dat->error_prob);
		free(dat);
	}
} /* free_data */

/**
 * Once data is loaded, build hash and optionally compress data.
 *
 * @param dat	data pointer
 * @param opt	options pointer
 * @return	error status
 */
int sync_data(data *dat, options *opt)
{
	UNUSED(opt);
	int err = NO_ERROR;
	
	/* rebuild a new hash table if previous hash table exist */
	if(dat->seq_count)
		delete_all(dat->seq_count);
	dat->seq_count = NULL;

	if ((err = build_hash(dat)))
		return err;

	return err;
} /* sync_data */

/* fill data object with provided dmat and qmat */
int fill_data(data *dat, data_t **dmat, data_t **qmat, unsigned int rlen, 
			size_t sample_size, unsigned int n_quality){

	int err = NO_ERROR;

	 dat->sample_size = sample_size;
     dat->max_read_length = rlen;
	 dat->max_read_position = rlen;
     dat->n_quality = n_quality;
     dat->dmat = dmat;
     dat->qmat = qmat;

    dat->error_prob = malloc(dat->n_quality * sizeof * dat->error_prob);
	if(!dat->error_prob)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,"dat.error_prob");

	for (unsigned int q = 0; q < dat->n_quality; q++)
		dat->error_prob[q] = error_prob(dat->fdata, q);

	dat->lengths = malloc(dat->sample_size * sizeof *dat->lengths);

	if (!dat->lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "data.lengths");

    for (size_t i = 0; i < dat->sample_size; ++i)
        dat->lengths[i] = dat->max_read_length;

	/* allocate the index array of reads */
	dat->read_idx = malloc(dat->fdata->n_reads * sizeof *dat->read_idx);

	if (!dat->read_idx)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "dat.read_idx");

	for (size_t i = 0; i < dat->fdata->n_reads; ++i)
		dat->read_idx[i] = i;

    if ((err = build_hash(dat)))
		return err;

	return NO_ERROR;
}/* fill_data */