/**
 * @file partition.c
 * @author Xiyu Peng
 *
 * Functions related to given partitions
 */

#include <ctype.h>

#include "io.h"
#include "error.h"
#include "constants.h"
#include "util.h"
#include "amplici.h"
#include "partition.h"
#include "io.h"
#include "libamplici.h"


/**
 * Run ampliCI cluster with an existing partition of the reads.
 *
 * @param opt	options object
 * @param dat	data object
 * @param mod	model object
 * @param ini	initialization object
 * @param ri	run information object
 * @error	error status
 */
int ampliCI_wpartition(options *opt, data *dat, model *mod, initializer *ini, run_info *ri)
{


	int err = NO_ERROR;
	int fxn_debug = opt->info;

	unsigned int K_space = opt->K;
	unsigned int ini_K = opt->K;
	unsigned int max_size = 0;

	//double abun_cutoff = opt->low_bound;

	/* read the partition file, start from 0 */
	/* need avoid NA in the partition file */
	if ((err = read_partition_file(opt->partition_file, ini->cluster_id,
							dat->sample_size)))
		return err;
	
	unsigned int max = 0;
	for (unsigned int i = 0; i < dat->sample_size; i++) 
		if (ini->cluster_id[i] > max)
			max = ini->cluster_id[i];

	unsigned int num_partition = max + 1;  // number of partitions 

	// realloc space for the number of partitions 
	unsigned int *cluster_size = realloc(ini->cluster_size,
						num_partition * sizeof *ini->cluster_size);
	
	if (!cluster_size)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"initializer::cluster_size");

	ini->cluster_size = cluster_size;

	for(unsigned int p = 0; p < num_partition; p++) 
		ini->cluster_size[p] = 0;

	for(unsigned int i = 0; i < dat->sample_size; i++)
		ini->cluster_size[ini->cluster_id[i]]++;

	for(unsigned int p = 0; p < num_partition; p++) 
		if (ini->cluster_size[p] > max_size)
			max_size = ini->cluster_size[p];

	/* alloc space for DNA sequence */
	data_t **dmat = NULL;
	data_t **qmat = NULL;

	dmat = malloc(max_size * sizeof *dmat);
	qmat = malloc(max_size *sizeof *qmat);

	if (!dmat || !qmat)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "partition: dmat & qmat");

	unsigned int sum_K = 0;
	// hash * hash_list = NULL;  // avoid duplicate 

	for (unsigned int p = 0; p < num_partition; p++){

	   //mmessage(INFO_MSG, NO_ERROR, "start the %d th subset\n",p);

		/* find the maximum size and reduce allocation */
		unsigned int size = ini->cluster_size[p];  // size of the subset

		/* get a subset of the data */

		/*
		data_t* rsub = NULL;
		data_t* qsub = NULL;
		
		rsub = malloc(dat->max_read_length * size * sizeof (*rsub));
		qsub = malloc(dat->max_read_length * size * sizeof (*rsub));

		if(!rsub ||!qsub)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "partition: rsub & qsub");
		
		*/
		unsigned int ss = 0;
		//unsigned char *rptr = rsub;
		//unsigned char *qptr = qsub;
		for(unsigned int i = 0; i < dat->sample_size; i++){
			if(ini->cluster_id[i] == p){
			   //memcpy(rptr, dat->dmat[i], dat->max_read_length * sizeof **dat->dmat);
			   //memcpy(qptr, dat->qmat[i], dat->max_read_length * sizeof **dat->qmat);
				dmat[ss] = dat->dmat[i];
				qmat[ss] = dat->qmat[i];
				ss++;
			   //rptr += dat->max_read_length;
			   //qptr += dat->max_read_length;
			}
		}

		//if(ss > size)
		//   return mmessage(ERROR_MSG, INTERNAL_ERROR, "error !");

			
		 /* prepare output for amplici_core */
		unsigned int K = 0;
		unsigned int *cluster_id = NULL;  // not needed here
		unsigned int *cluster_size = NULL;  // not needed here
		unsigned char *seeds = NULL;
		unsigned int *seeds_length = NULL;
		double *ll = NULL;
		double *abun = NULL;

		// mmessage(INFO_MSG, NO_ERROR, "stat amplicon clustering \n",p);

		err = amplici_core(dmat, qmat, (size_t) size, opt->low_bound,dat->max_read_length, 
				(char*) opt->error_profile_name, dat->fdata->max_quality, dat->fdata->min_quality, 
				&seeds, &seeds_length, &cluster_id, 
				&cluster_size, &K, &ll, &abun);

		// OK when K = 0
		if(err && K > 0)
			return err;

		//mmessage(INFO_MSG, NO_ERROR, "finish amplicon clustering \n",p);


	   // mmessage(INFO_MSG, NO_ERROR, "start store seeds\n",p);

		/* store seeds */
		//int first;
		for (unsigned int k = 0; k < K; ++k) {
			//mmessage(INFO_MSG, NO_ERROR, "store the %d th hap\n",k);
			//mmessage(INFO_MSG, NO_ERROR, "seeds_length =  %d \n",seeds_length[k]);
			//for (size_t j = 0; j < seeds_length[k]; ++j) {
			//		fprintf(stderr, "%c", xy_to_char[(int) seeds[k*dat->max_read_length+j]]);
			//}
			hash *new;
			HASH_FIND(hh, dat->seq_count, &seeds[k*dat->max_read_length], seeds_length[k] * sizeof (unsigned char), new);
			if(new->seeds == 0){   // avoid duplicate
				new->seeds = 1;
				sum_K ++;
			}
			//if(!new)
			//	mmessage(INFO_MSG, NO_ERROR, "find it %d !\n",new->count);

			/* 
			first = 0;
			if(abun[k] >= abun_cutoff)
				first = add_sequence(&hash_list, &seeds[k*dat->max_read_length],
					seeds_length[k], k, &err);
			mmessage(INFO_MSG, NO_ERROR, "fisrt = %d \n",first);
			if(err) return err;
			if(first){
				mmessage(INFO_MSG, NO_ERROR, "store the %d th hap\n",sum_K);
				memcpy(ini->seeds[sum_K], &seeds[k*dat->max_read_length], seeds_length[k] * sizeof *seeds);
				ini->seed_lengths[sum_K] = seeds_length[k];
				sum_K ++;
			}
			*/
		}

	   //mmessage(INFO_MSG, NO_ERROR, "finish store seeds\n",p);

		
		if(cluster_id)free(cluster_id);
		if(cluster_size) free(cluster_size);
		if(seeds) free(seeds);
		if(seeds_length) free(seeds_length);
		if(ll) free(ll);
		if(abun) free(abun);
		//if(rsub) free(rsub);
		//if(qsub) free(qsub);
	}

	if(dmat)free(dmat);
	if(qmat)free(qmat);

	 /* realloc the space for seeds for new K */
	if (sum_K  > K_space) {
			debug_msg(DEBUG_III, fxn_debug, "begin reallocation\n");
				 
			/* for seeds */
			if((err = realloc_seeds(ini, dat->max_read_length, K_space, sum_K)))
				return err;
	
			debug_msg(DEBUG_III, fxn_debug,
						"Finish reallocation\n");
	}
	opt->K = sum_K;

	/* construct seeds set */
	hash *s;
	unsigned int cnt = 0;
	for (s = dat->seq_count; s != NULL; s = s->hh.next) {
			if(s->seeds == 1){
				memcpy(ini->seeds[cnt], s->sequence,
					dat->max_read_length * sizeof **ini->seeds);
				ini->seed_lengths[cnt] = dat->max_read_length;
				cnt++;
			}
	}
	if(sum_K != cnt)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "sum_K != cnt");

	/* output */
	char *outfile_hap = NULL;

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
		
		fprint_haplotypes_abun(fp2,&ini->seeds[0], opt->K,
			&dat->max_read_length, opt->p_threshold, "H", 
				NULL, NULL, NULL);

		fclose(fp2);

		mmessage(INFO_MSG, NO_ERROR, "Output the final haplotype fasta "
					"file: %s \n", opt->outfile_fasta);

	return err;
}/* ampliCI_wpartition */


/** 
 * Read the partition file.
 *
 * @param filename	partition file
 * @param cluster_id	array for storing id
 * @param sample_size	length of the array  
 * @return		error status
 */
int read_partition_file(char const * const filename, unsigned int *cluster_id, unsigned int sample_size)
{
	int err = NO_ERROR;

	FILE *fp = fopen(filename, "r");

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	/* detect if AmpliCI out file */
	char c = fgetc(fp);
	ungetc(c, fp);

	//mmessage(DEBUG_MSG, NO_ERROR, "Read char '%c'.\n", c);

	if (!isdigit(c)) {	/* look for assignments: */
		char word[12];

		do {
			next_word(fp, word, 11, ':');
			//mmessage(DEBUG_MSG, NO_ERROR, "Word '%s'\n", word);
			if (feof(fp))
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
								filename);
			if (!strcmp(word, "assignments")) {
				fgetc(fp);	/* : */
				break;
			} else {
				fforward(fp, c, '\n');
			}
		} while (1);
	}

	err = fread_uints(fp, cluster_id, sample_size);
	
	fclose(fp);
	return err;
}/* read_partition_file */
