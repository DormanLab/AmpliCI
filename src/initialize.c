/**
 * @file initialize.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 * 
 * Initialize ampliclust.
 *
 * TODO
 * - initialization is unaware of \par options.background_model
 * - initialization clusters haplotypes of length \par data.min_read_length;
 *   initializes remaining positions in a foolish way
 *
 * DONE
 * X initialize from set of seeds
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <string.h>
#include <stdlib.h>

#include "initialize.h"
#include "data.h"
#include "options.h"

/**
 * Create initializer object.
 *
 * @param ini	initializer object
 * @param dat	data object
 * @param opt	options object
 * @return	error status
 */
int make_initializer(initializer **ini, data *dat, options *opt,fastq_data *fqdf)
{
	int err = NO_ERROR;
	initializer *in;
	*ini = malloc(sizeof **ini);

	if (*ini == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer object");
	in = *ini;

	in->K = opt->K;
	in->cluster_id = NULL;

	in->seed_idx = calloc(opt->K, sizeof *in->seed_idx);

	if (!in->seed_idx)
		return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
			MEMORY_ALLOCATION, "initializer::seed_idx");

	in->criterion = malloc(opt->K * sizeof *in->criterion);

	if (!in->criterion)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer::criterion");
	in->cluster_size = malloc((opt->K)
		* sizeof *in->cluster_size);
	
	if (!in->cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer::cluster_size");

	in->optimal_total = INFINITY;

	/* allocate room for k-modes initialization */
	in->seeds = malloc(opt->K * sizeof *in->seeds);
	in->seed_lengths = calloc(opt->K, sizeof *in->seed_lengths);

	if (!in->seeds || !in->seed_lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "initializer.seeds");

	/* [TODO] allocate seeds in one block: easier to free and realloc:
	 * data_t *dptr = malloc(dat->max_read_length * opt->K * sizeof **in_seeds);
	 * for (size_t i = 1; i < opt->K; ++i) {
	 *	in->seeds[i] = dptr;
	 *	dptr += dat->max_read_length;
	 * }
	 * //... use in->seeds ... //
	 * if (in->seeds) {
	 *	if (in->seeds[0])
	 *		 free(in->seeds[0])
	 *	free(in->seeds);
	 * }
	 * */
	for (size_t i = 0; i < opt->K; i++) {
		in->seeds[i] = calloc(dat->max_read_length, sizeof **in->seeds);
		if (in->seeds[i] == NULL )
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"initializer.seeds");
	}

	/* allocate room for cluster assignments */
	in->cluster_id = calloc(dat->sample_size, sizeof *in->cluster_id);

	if (!in->cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer.cluster_id");

	if(fqdf){
		for (unsigned int k = 0; k < opt->K; ++k){
			memcpy(in->seeds[k],
				&fqdf->reads[k * fqdf->n_max_length],
				fqdf->n_max_length * sizeof *fqdf->reads);
			in->seed_lengths[k] = fqdf->n_max_length;
		//	mmessage(INFO_MSG, NO_ERROR, "%u\n", k);
		//	mmessage(INFO_MSG, NO_ERROR, "Read %2u: %.*s\n", k,
		//		in->seed_lengths[k], display_sequence(in->seeds[k],
		//		in->seed_lengths[k], XY_ENCODING));
			//display_sequence(in->seeds[k], in->seed_lengths[k], XY_ENCODING);			
		}
		//mmessage(INFO_MSG, NO_ERROR, "number of haplotypes is "
		//	"%u!\n", fqdf->n_reads);
	}

	in->abun_true = NULL;
	in->p = NULL;
	in->H = NULL;
	in->H_abun = NULL;
	in->H_ee = NULL;
	in->H_pvalue = NULL;
	in->e_trans = NULL;
	in->self_trans = NULL;
	in->nw_mismatch = NULL;
	in->nw_indels = NULL;

	in->err_cnt = NULL;

	in->uniq_seq_idx = NULL;
	in->uniq_seq_count = NULL;
	in->reads_uniq_id = NULL;

	if (!dat->fdata->empty)
		err = sync_initializer(in, dat);

	return err;
} /* make_initializer */


/**
 * Finish initializer setup once data available.
 *
 * @param ini	initializer pointer
 * @param dat	data pointer
 * @return	error status
 */
int sync_initializer(initializer *ini, data *dat)
{
	/* allocate room for index and count array of unique sequences */
	ini->uniq_seq_idx = calloc(dat->hash_length, sizeof *ini->uniq_seq_idx);
	ini->uniq_seq_count = calloc(dat->hash_length,
						sizeof *ini->uniq_seq_count);
	ini->reads_uniq_id = calloc(dat->sample_size,sizeof *ini->reads_uniq_id);

	if (!ini->uniq_seq_idx || !ini->uniq_seq_count || !ini->reads_uniq_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"initializer.uniq_seq");

	/* hash table should be sorted in an order of decreasing */
	if (store_index(dat->seq_count, ini->reads_uniq_id, ini->uniq_seq_idx, 
		dat->hash_length,dat->sample_size)
		|| store_count(dat->seq_count, ini->uniq_seq_count,
							dat->hash_length))
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"store_index() or store_count()");

	return NO_ERROR;
} /* sync_initializer */


/**
 * Reallocate initializer struct for a different K or sample_size.  It is
 * assumed that options::K and data::sample_size contain the new numbers.
 *
 * @param ini		pointer to initialization object
 * @param dat		pointer to data object
 * @param opt		pointer to options object
 * @return          error status
 */
int realloc_initializer(initializer *ini, data *dat, options *opt)
{
	/* reallocate the room for a new K*/
	/* seed_idx */
	/* claculate the previous K */
	unsigned int pre_K = ini->K;

	if(pre_K <2)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,"realloc.initializer");
	ini->K = opt->K;

	size_t *seed_idx = realloc(ini->seed_idx,
					opt->K * sizeof *ini->seed_idx);

	if (!seed_idx )
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.seed_idx");

	ini->seed_idx = seed_idx;

	/* criterion */
	double *criterion = realloc(ini->criterion,
					opt->K * sizeof *ini->criterion);

	if (!criterion)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.criterion");

	ini->criterion = criterion;

	/* seeds, best_modes,seed_lengths*/
	/* free previous points */
	/* [TODO] easier if you allocate ini->{seeds, best_modes, etc} in one block (see other TODOs) */

	if (ini->seeds) {
		for (unsigned int i = 0; i < pre_K; ++i){
			if (ini->seeds[i])
				free(ini->seeds[i]);
		}
		free(ini->seeds);

	}
	ini->seeds = NULL;

	data_t **seeds = malloc(opt->K * sizeof *ini->seeds);

	unsigned int *seed_lengths = realloc(ini->seed_lengths,
				opt->K * sizeof *ini->seed_lengths);

	if (!seeds || !seed_lengths)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"reallloc.initializer.seeds");

	ini->seeds = seeds;
	ini->seed_lengths = seed_lengths;

	for (unsigned int i = 0; i < opt->K; ++i) {
		if (!dat->max_read_length)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds");
		ini->seeds[i] = calloc(dat->max_read_length, sizeof *ini->seeds);
		if (ini->seeds[i] == NULL )
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds");
	}

	/* cluster_size */
	unsigned int *cluster_size = realloc(ini->cluster_size,
		(opt->K) * sizeof *ini->cluster_size);
	if (!cluster_size )
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.cluster_size");

	ini->cluster_size = cluster_size;
	
	/* cluster_id */
	unsigned int *cluster_id = realloc(ini->cluster_id,
		dat->sample_size * sizeof *ini->cluster_id);

	if (!cluster_id )
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.initializer.cluster_id");

	ini->cluster_id = cluster_id;

	/* index and count array of unique sequences */
	size_t *uniq_seq_idx = realloc(ini->uniq_seq_idx,
		dat->hash_length * sizeof *ini->uniq_seq_idx);
	unsigned int *uniq_seq_count = realloc(ini->uniq_seq_count,
		dat->hash_length * sizeof *ini->uniq_seq_count);
	unsigned int *reads_uniq_id = realloc(ini->reads_uniq_id,
		dat->sample_size * sizeof *ini->reads_uniq_id);

	if (!uniq_seq_idx || !uniq_seq_count || !reads_uniq_id)
		return mmessage(ERROR_MSG,MEMORY_ALLOCATION,"realloc.initializer.uniq_seq");

	ini->uniq_seq_idx = uniq_seq_idx;
	ini->uniq_seq_count = uniq_seq_count;
	ini->reads_uniq_id = reads_uniq_id;

	if(ini->abun_true) free(ini->abun_true);
	if(ini->p) free(ini->p);
	if(ini->H) free(ini->H);
	if(ini->H_abun) free(ini->H_abun);
	if(ini->e_trans) free(ini->e_trans);
	if(ini->self_trans) free(ini->self_trans);
	if (ini->nw_mismatch) free(ini->nw_mismatch);
	if (ini->nw_indels) free(ini->nw_indels);

	ini->abun_true = NULL;
	ini->p = NULL;
	ini->H = NULL;
	ini->H_abun = NULL;
	ini->e_trans = NULL;
	ini->self_trans = NULL;
	ini->nw_mismatch = NULL;
	ini->nw_indels = NULL;

	/* hash table should have been sorted in an order of decreasing */
	if (store_index(dat->seq_count, ini->reads_uniq_id,ini->uniq_seq_idx, 
		dat->hash_length,dat->sample_size)
		|| store_count(dat->seq_count, ini->uniq_seq_count,
							dat->hash_length))
		return mmessage(ERROR_MSG, INTERNAL_ERROR,
					"store_index() or store_count()");

	ini->optimal_total = INFINITY;
	
	return NO_ERROR;
}/* realloc_initializer */


/**
 * Read initialization information.  Initialization file can be a fasta file
 * containing haplotypes, 
 *
 *
 * @param filename	name of file containing initialization information
 * @param dat		data object pointer
 * @param opt		options object pointer
 * @param ini		initializer object pointer
 * @return		error status
 */
int read_initialization_file(char const * const filename, data *dat,
	options *opt, fastq_data **fqdf){

	int fxn_debug = opt->info;	//DEBUG_III;	//
	int err = NO_ERROR;
	FILE *fp = fopen(filename, "r");

	debug_msg(DEBUG_III, fxn_debug, "Opening file '%s'\n", filename);

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	char c = fgetc(fp);

	rewind(fp);

	//debug_msg(DEBUG_III, fxn_debug, "First character '%c'\n", c);

	
	debug_msg(DEBUG_III, fxn_debug, "entering fasta read\n");
	//mmessage(INFO_MSG, NO_ERROR, "Read the haplotype set: "
	//			"%s\n", opt->initialization_file);

	/* read in fasta-formatted haplotypes */
	fastq_data *fqd = NULL;
	
	fastq_options fop = {.read_encoding = XY_ENCODING};
	if ((err = fread_fastq(fp, fqdf, &fop))) {
		debug_msg(DEBUG_III, fxn_debug, "err=%d\n", err);
		fclose(fp);
		return err;
	}

	fqd = *fqdf;

	mmessage(INFO_MSG, NO_ERROR, "finished reading %u sequences "
		"in %s format\n", fqd->n_reads,
		fqd->file_type == FASTA_FILE ? "fasta" : "fastq");

	if (fqd->n_lengths || fqd->n_max_length != dat->max_read_length) {
		fclose(fp);
		mmessage(INFO_MSG, NO_ERROR, "haplotypes in '%s' must be same "
							"length.\n", filename);
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
						"invalid input haplotype set");
	}

	//mmessage(INFO_MSG, NO_ERROR, "number of haplotypes in '%s' is "
	//		"%u!\n", filename, fqd->n_reads);

	/* KSD: Why not just set opt->K = fqd->n_reads? You just need to call this function earlier, like during make_initializer(). */
	/*
	if (fqd->n_reads != opt->K) {	
		fclose(fp);
		mmessage(INFO_MSG, NO_ERROR, "number of haplotypes in '%s' is "
			"%u! Rerun with -k %u !\n", filename, fqd->n_reads,
								fqd->n_reads);
		return mmessage(ERROR_MSG, INVALID_USER_INPUT, "invalid input K");
	}
	*/

	/* move to make_initializer */
	/*
	for (unsigned int k = 0; k < opt->K; ++k){
		memcpy(ini->seeds[k],
			&fqd->reads[k * fqd->n_max_length],
			fqd->n_max_length * sizeof *fqd->reads);
		ini->seed_lengths[k] = fqd->n_max_length;			
	}
	*/

	//free_fastq(fqd);

	fclose(fp);
	return err;
} /* read_initialization_file */

void free_initializer(initializer *ini, options *opt)
{
	if (ini) {
		if (ini->cluster_id) free(ini->cluster_id);
		if (ini->seed_idx) free(ini->seed_idx);

		if (ini->criterion) free(ini->criterion);

		if (ini->cluster_size) free(ini->cluster_size);

		unsigned int size;
		if (ini->seeds) {
			size = opt->run_amplici?opt->K_space:opt->K;
			for (size_t i = 0; i < size; ++i)
			{
				if(ini->seeds[i])
					free(ini->seeds[i]);
			}
			free(ini->seeds);
		}

		if (ini->seed_lengths) free(ini->seed_lengths);
		if (ini->uniq_seq_count) free(ini->uniq_seq_count);
		if (ini->uniq_seq_idx) free(ini->uniq_seq_idx);
		if (ini->reads_uniq_id) free(ini->reads_uniq_id);
		if (ini->abun_true) free(ini->abun_true);
		if (ini->p) free(ini->p);
		if (ini->H) free(ini->H);
		if (ini->H_abun) free(ini->H_abun);
		if (ini->H_ee) free(ini->H_ee);
		if (ini->H_pvalue) free(ini->H_pvalue);
		if(ini->e_trans) free(ini->e_trans);
		if (ini->self_trans) free (ini->self_trans);
		if (ini->nw_mismatch) free(ini->nw_mismatch);
		if (ini->nw_indels) free(ini->nw_indels);
		if (ini->err_cnt) free(ini->err_cnt); 

		free(ini);
	}
} /* free_initializer */
