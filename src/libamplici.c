/**
 * @file libamplici.c
 * @author Xiyu Peng
 *
 * A C library for clustering amplicon sequences
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "ampliclust.h"
#include "data.h"
#include "model.h"
#include "options.h"
#include "initialize.h"
#include "fastq.h"
#include "io.h"
#include "amplici.h"
#include "error_est.h"
#include "libamplici.h"

/**
 * Cluster amplicon sequences with a given fastq file
 * 
 * @param fastq_file	    the name of the fastq file    [REQUIRED]
 * @param error_profile   the name of the error profile [OPTION]
 * @param seeds           pointer to seeds (haplotypes) (K x l) [OUTPUT]
 * @param seeds_length    pointer to seeds length (K x 1) [OUTPUT]
 * @param cluster_id      pointer to reads cluster index (N x 1) [OUTPUT]
 * @param cluster_size    pointer to cluster size (K x 1) [OUTPUT]
 * @param K               pointer to number of clusters  [OUTPUT]
 * @param sample_size     pointer to sample size (N)  [OUTPUT]
 * @param max_read_length poonter to maximum read length (l) [OUTPUT]
 *              (The kth (in [0,1,2,...K-1]) haplotype starts at seeds[k*l])
 * @return err status
 **/
int amplici_wfile(char *fastq_file, char *error_profile_name, unsigned char **seeds, unsigned int **seeds_length, 
               unsigned int **cluster_id, unsigned int **cluster_size, unsigned int *K, size_t *sample_size,
               unsigned int* max_read_length){

  int err = NO_ERROR;		/* error state */
	options *opt = NULL;		/* run options */
	data *dat = NULL;		/* data object */
	model *mod = NULL;		/* model object */
	fastq_options *fqo = NULL;	/* fastq file options */
	initializer *ini = NULL;	/* initializer */

  if (!fastq_file)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"No input data");

	/* parse command line */
	if ((err = make_options(&opt)))
		goto AMPLICI_CLEAR2;

	opt->K = opt->K_space; // set intial max K
  opt->fastq_file = fastq_file;
  opt->error_profile_name = error_profile_name;

	/* make data object */
	if ((err = make_data(&dat, opt)))
		goto AMPLICI_CLEAR2;

	if ((err = make_fastq_options(&fqo)))
		goto AMPLICI_CLEAR2;

	/* encode nucleotides in 2-bits: error raised if ambiguous bases */
	fqo->read_encoding = XY_ENCODING;

	/* read sequence data */
	if (opt->fastq_file && (err = read_fastq(opt->fastq_file,
						&dat->fdata, fqo)))
		goto AMPLICI_CLEAR2;

	/* with data now loaded, can polish off data object */
	if ((err = sync_state(dat, opt)))
		goto AMPLICI_CLEAR2;

	/* create model
	 * [TODO] reasonable defaults, but uses data for binned quality models
	 */
	if ((err = make_model(&mod, dat, opt)))
		goto AMPLICI_CLEAR2; 

  //mod->error_profile = error_profile;

	/* make initializer */
	if ((err = make_initializer(&ini, dat, opt, NULL, NULL)))
		goto AMPLICI_CLEAR2;

  if (error_profile_name)
    opt->use_error_profile = 1;

  if ((err = haplotype_selection(opt, dat, mod, ini, opt->K_max)))
    return err;

  /* assgin reads */
  assign_clusters(mod->eik, opt->K, dat->sample_size, ini->cluster_size,
                  ini->cluster_id, 1);

  /* copy to seeds, cluster_id, cluster_size */
  *K = opt->K;
  *sample_size = dat->sample_size;
  *cluster_id = ini->cluster_id;
  *cluster_size = ini->cluster_size;
  *seeds_length = ini->seed_lengths;
  *seeds = ini->seeds[0];  // check carefully 
  *max_read_length = dat->max_read_length;

  /* avoid to be freed */
  /* input */
  mod->error_profile = NULL;
  opt->fastq_file = NULL;

  /* output */
  ini->cluster_id = NULL;
  ini->cluster_size = NULL;
  ini->seeds[0] = NULL;
  ini->seed_lengths = NULL;
  

  AMPLICI_CLEAR2:
  if (dat)
    free_data(dat);
  if (ini)
    free_initializer(ini, opt);
  if (mod)
    free_model(mod);
  if (opt)
    free_options(opt);
  if (fqo)
		free_fastq_options(fqo);
  
  return NO_ERROR;
}/*  amplici_wfile  */


/* Another function to call haplotype selection [UNDER DEVELOPMENT] */
int amplici_core(data_t **dmat, data_t **qmat, size_t sample_size, unsigned int rlen, double *error_profile, unsigned int n_quality, 
            data_t ***seeds, unsigned int **seeds_length, unsigned int **cluster_id, unsigned int **cluster_size, unsigned int *K)
{

  int err = NO_ERROR;

  options *opt = NULL;     /* run options */
  data *dat = NULL;        /* data object */
  model *mod = NULL;       /* model object */
  initializer *ini = NULL; /* initializer */

  /* make option object */
  if ((err = make_options(&opt)))
    goto AMPLICI_CLEAR;

  opt->K = opt->K_space; // set intial max K

  /* make data object */
  if ((err = make_data(&dat, opt)))
    goto AMPLICI_CLEAR;

  /* fill missing information in data object */
  if ((err = fill_data(dat, dmat, qmat, rlen, sample_size, n_quality)))
    goto AMPLICI_CLEAR;

  if ((err = make_model(&mod, dat, opt)))
    goto AMPLICI_CLEAR;
  mod->error_profile = error_profile;

  if ((err = make_initializer(&ini, dat, opt, NULL, NULL)))
    goto AMPLICI_CLEAR;

  if (error_profile)
    opt->use_error_profile = 1;

  if ((err = haplotype_selection(opt, dat, mod, ini, opt->K_max)))
    return err;

  /* assgin reads */
  assign_clusters(mod->eik, opt->K, dat->sample_size, ini->cluster_size,
                  ini->cluster_id, 1);

  /* copy to seeds, cluster_id, cluster_size */
  *K = opt->K;
  *cluster_id = ini->cluster_id;
  *cluster_size = ini->cluster_size;
  *seeds = ini->seeds;
  *seeds_length = ini->seed_lengths;

  /* avoid to be freed */
  dat->dmat = NULL;
  dat->qmat = NULL;
  mod->error_profile = NULL;
  ini->cluster_id = NULL;
  ini->cluster_size = NULL;
  ini->seeds = NULL;
  ini->seed_lengths = NULL;

AMPLICI_CLEAR:
  if (dat)
    free_data(dat);
  if (ini)
    free_initializer(ini, opt);
  if (mod)
    free_model(mod);
  if (opt)
    free_options(opt);

  return err;
}/* amplici_core */