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
int make_initializer(initializer **ini, data *dat, options *opt,fastq_data *fqdf,fastq_data *fqdfu)
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

	in->seed_idx = NULL;  //calloc(opt->K, sizeof *in->seed_idx);

	//if (!in->seed_idx)
	//	return message(stderr, __FILE__, __func__, __LINE__, ERROR_MSG,
	//		MEMORY_ALLOCATION, "initializer::seed_idx");

	in->criterion = NULL; // malloc(opt->K * sizeof *in->criterion);

	//if (!in->criterion)
	//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
	//		"initializer::criterion");

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

	data_t *dptr = malloc(dat->max_read_length * opt->K * sizeof **in->seeds);
	if(!dptr)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,"initializer.seeds");
	for (size_t i = 0; i < opt->K; i++) {
		in->seeds[i] = dptr;
		dptr += dat->max_read_length;
	}

	/* allocate room for cluster assignments */
	in->cluster_id = calloc(dat->sample_size, sizeof *in->cluster_id);

	if (!in->cluster_id)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer.cluster_id");

	if(fqdf){
		if (fqdf->n_lengths || fqdf->n_max_length != dat->max_read_length) 
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
						"invalid input haplotype set. The length of input haplotypes in '%s' must be "
							"%i.\n",opt->initialization_file, dat->max_read_length);
		
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

	if (dat->dmat && dat->qmat)
		err = sync_initializer(in, dat);

	/* UMI information */
	in->seeds_UMI = NULL;
	in->seeds_hash = NULL;
	in->UMIs_hash = NULL;
	in->reads_hap_id = NULL;
	in->reads_umi_id = NULL;

	int fxn_debug = ABSOLUTE_SILENCE;
 
	if (opt->UMI_length)
	{	/* seeds_UMI */
		in->seeds_UMI = malloc(opt->UMI_length * opt->K_UMI * sizeof *in->seeds_UMI);
		if (!in->seeds_UMI)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"initializer.seeds_UMI");
		if (fqdfu)
		{
			if (fqdfu->n_lengths || fqdfu->n_max_length != opt->UMI_length)
				return mmessage(ERROR_MSG, INVALID_USER_INPUT,
								"invalid input UMI set. UMIs in '%s' must be "
								"%i.\n",opt->initialization_UMI, opt->UMI_length);

			memcpy(in->seeds_UMI, fqdfu->reads,
				   opt->K_UMI * fqdfu->n_max_length * sizeof *fqdfu->reads);
		}

		/* seeds_hash */
		for (size_t k = 0; k < opt->K; ++k) {
			int first = add_sequence(&in->seeds_hash, in->seeds[k],
							in->seed_lengths[k], k,&err);
			if(!first)
				mmessage(WARNING_MSG, INVALID_USER_INPUT,
								"Duplicate sequences in '%s'.\n",
								opt->initialization_file);
			if(err) return err;
		}

		/* reads_hap_id (error free) */
		in->reads_hap_id = malloc(dat->sample_size * sizeof *in->reads_hap_id);
		if(!in->reads_hap_id)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,"initializer.reads_hap_id");

		for (size_t i = 0; i < dat->sample_size; ++i){
			hash *new;
			HASH_FIND( hh, in->seeds_hash, dat->dmat[i], dat->lengths[i] * sizeof *dat->dmat[i], new);
			if(new){ in->reads_hap_id[i] = new->idx; }else {in->reads_hap_id[i] = -1;}
			debug_msg(DEBUG_I, fxn_debug, "%5d,",in->reads_hap_id[i]);
		}
		debug_msg(DEBUG_I, fxn_debug, "\n");

		/* UMIs_hash */
		for (size_t k = 0; k < opt->K_UMI; ++k) {
			int first = add_sequence(&in->UMIs_hash, &in->seeds_UMI[k*opt->UMI_length],
							opt->UMI_length, k, &err);
			if(!first)
				mmessage(WARNING_MSG, INVALID_USER_INPUT,
								"Duplicate sequences in '%s'.\n",
								opt->initialization_UMI);
			if(err) return err;
		}

		/* reads_umi_id (error free) */
		in->reads_umi_id = malloc(dat->sample_size * sizeof *in->reads_umi_id);
		if(!in->reads_umi_id)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,"initializer.reads_hap_id");

		for (size_t i = 0; i < dat->sample_size; ++i){
			hash *newU;
			HASH_FIND(hh, in->UMIs_hash, dat->dmatU[i], opt->UMI_length * sizeof *dat->dmatU[i], newU);
			if(newU){ in->reads_umi_id[i] = newU->idx; }else {in->reads_umi_id[i] = -1;}
			debug_msg(DEBUG_I, fxn_debug, "%5d,",in->reads_umi_id[i]);
		}
		debug_msg(DEBUG_I, fxn_debug, "\n");


		/* H_abun */
		in->H_abun = calloc(opt->K, sizeof (*in->H_abun));
		if(!in->H_abun)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "initializer.H_abun");	
	}


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
	unsigned int preK = ini->K;
	unsigned int err;

	if(preK <2)
		return mmessage(ERROR_MSG, INTERNAL_ERROR,"realloc.initializer");
	ini->K = opt->K;
	size_t *seed_idx = NULL;
	double *criterion = NULL;

	//size_t *seed_idx = realloc(ini->seed_idx,
	//				opt->K * sizeof *ini->seed_idx);

	//if (!seed_idx )
	//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
	//		"realloc.initializer.seed_idx");

	//ini->seed_idx = seed_idx;

	/* criterion */
	//double *criterion = realloc(ini->criterion,
	//				opt->K * sizeof *ini->criterion);

	//if (!criterion)
	//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
	//		"realloc.initializer.criterion");

	//ini->criterion = criterion;

	/* seeds, best_modes,seed_lengths*/
	/* free previous points */
	/* [TODO] easier if you allocate ini->{seeds, best_modes, etc} in one block (see other TODOs) */
	if((err = realloc_seeds(ini, dat->max_read_length, preK, opt->K)))
		return err;
	/* 
	data_t **seeds = realloc(ini->seeds, opt->K * sizeof *ini->seeds);

	unsigned int *seed_lengths = realloc(ini->seed_lengths,
				opt->K * sizeof *ini->seed_lengths);

	if (!seeds || !seed_lengths){
		if (seeds) free(seeds);
		if (seed_lengths) free(seed_lengths);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"reallloc.initializer.seeds");
	}

	ini->seeds = seeds;
	ini->seed_lengths = seed_lengths;

	data_t *dptr = realloc(ini->seeds[0], dat->max_read_length * opt->K * sizeof **ini->seeds);
	if (!dptr)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"reallloc.initializer.seeds");
	for (size_t i = 0; i < opt->K; i++) {
		ini->seeds[i] = dptr;
		dptr += dat->max_read_length;
	} */

	/*
	for (unsigned int i = 0; i < opt->K; ++i) {
		if (!dat->max_read_length)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds");
		ini->seeds[i] = calloc(dat->max_read_length, sizeof *ini->seeds);
		if (ini->seeds[i] == NULL )
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.initializer.seeds");
	}
	*/

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
	//if(ini->p) free(ini->p);
	if(ini->H) free(ini->H);
	if(ini->H_abun) free(ini->H_abun);
	if(ini->e_trans) free(ini->e_trans);
	if(ini->self_trans) free(ini->self_trans);
	if (ini->nw_mismatch) free(ini->nw_mismatch);
	if (ini->nw_indels) free(ini->nw_indels);

	ini->abun_true = NULL;
	// ini->p = NULL;
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
	
	return err;
}/* realloc_initializer */


/**
 * reallloc ini->seeds, ini->seeds_length for larger K 
 * 
 * @param ini  initializer object
 * @param max_read_length space for each seed
 * @param preK  Previous malloc spaces for K haplotypes
 * @param K     realloc space for K haplotypes
 * 
 * return 
 */
int realloc_seeds(initializer *ini, unsigned int max_read_length, unsigned int preK, unsigned int K){
	
	UNUSED(preK);
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
	// if (ini->seeds[0] == dptr)  s = preK; [Do not understand why it is a bug]

	for (size_t k = s; k < K; k++) {
		ini->seeds[k] = dptr;
		dptr += max_read_length;
	}

	return NO_ERROR;
}/* realloc_seeds */


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
int read_initialization_file(char const * const filename, fastq_data **fqdf, int info){

	int fxn_debug = info;	//DEBUG_III;	//
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

	/* remove to support haplotypes with multiple length */
	/* 
	if (fqd->n_lengths || fqd->n_max_length != dat->max_read_length) {
		fclose(fp);
		mmessage(INFO_MSG, NO_ERROR, "haplotypes in '%s' must be same "
							"length.\n", filename);
		return mmessage(ERROR_MSG, INVALID_USER_INPUT,
						"invalid input haplotype set");
	}
	*/

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

		if (ini->seeds) {
			if (ini->seeds[0])
				free(ini->seeds[0]);	
			free(ini->seeds);
			ini->seeds = NULL;
		}
		if (ini->seed_lengths) free(ini->seed_lengths);
		if (ini->uniq_seq_count) free(ini->uniq_seq_count);
		if (ini->uniq_seq_idx) free(ini->uniq_seq_idx);
		if (ini->reads_uniq_id) free(ini->reads_uniq_id);
		if (ini->abun_true) free(ini->abun_true);
		// if (ini->p) free(ini->p);
		if (ini->H) free(ini->H);
		if (ini->H_abun) free(ini->H_abun);
		if (ini->H_ee) free(ini->H_ee);
		if (ini->H_pvalue) free(ini->H_pvalue);
		if(ini->e_trans) free(ini->e_trans);
		if (ini->self_trans) free (ini->self_trans);
		if (ini->nw_mismatch) free(ini->nw_mismatch);
		if (ini->nw_indels) free(ini->nw_indels);
		if (ini->err_cnt) free(ini->err_cnt); 
		if(ini->seeds_UMI) free(ini->seeds_UMI);
		if (ini->seeds_hash) delete_all(&ini->seeds_hash);
		if (ini->UMIs_hash) delete_all(&ini->UMIs_hash);
		if (ini->reads_hap_id)free(ini->reads_hap_id);
		if (ini->reads_umi_id) free(ini->reads_umi_id);

		free(ini);
	}
} /* free_initializer */

