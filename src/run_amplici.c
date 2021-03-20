/**
 * @file run_amplici.c
 * @author Xiyu Peng
 * @author Karin S. Dorman
 *
 * Cluster amplicon sequences.
 *
 * TODO
 * - rewrite cluster_amplicons() so that writing a file of results is optional
 * - handle error indel in reads: very rare in Illumina
 * - unsigned char data_t (previous) VS uint8_t data_t (current);
 * - remove curses entirely?
 * - provide fqmorph for people to get rid of reads with ambiguous bases
 *
 * DONE
 * X created initializer object to store initialization information
 *   separate from data and model
 * X change default to NOT use curses
 * X split into separate files for options, model, data, initalization, etc.
 * X estimate error rates (& compare to estimated from true error data)
 * X allow different read lengths
 * X rewrite message handling to use global debug level
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
#include "amplici_umi.h"
#include "partition.h"

int main(int argc, const char **argv)
{
	int err = NO_ERROR;		/* error state */
	options *opt = NULL;		/* run options */
	data *dat = NULL;		/* data object */
	model *mod = NULL;		/* model object */
	fastq_options *fqo = NULL;	/* fastq file options */
	initializer *ini = NULL;	/* initializer */
	run_info *ri=NULL;		/*run_info object */
	fastq_data *fqdf = NULL;    /* initialize fasta file */
	fastq_data *fqdfu = NULL;    /* initialize fasta UMI file */

	/* parse command line */
	if ((err = make_options(&opt)))
		goto CLEAR_AND_EXIT;

	if ((err = parse_options(opt, argc, argv)))
		goto CLEAR_AND_EXIT;

	/* make data object */
	if ((err = make_data(&dat, opt)))
		goto CLEAR_AND_EXIT;

	if ((err = make_fastq_options(&fqo)))
		goto CLEAR_AND_EXIT;

	/* encode nucleotides in 2-bits: error raised if ambiguous bases */
	fqo->read_encoding = XY_ENCODING;

	/* read sequence data */
	if (opt->fastq_file && (err = read_fastq(opt->fastq_file,
						&dat->fdata, fqo)))
		goto CLEAR_AND_EXIT;

	/* with data now loaded, can polish off data object */
	if ((err = sync_state(dat, opt)))
		goto CLEAR_AND_EXIT;

	/* read initialization file  into fqdf */
	if (opt->initialization_file) {
		if ((err = read_initialization_file(opt->initialization_file,
							&fqdf, opt->info)))
			goto CLEAR_AND_EXIT;
		opt->K = fqdf->n_reads;
	}

	/* read initialization UMI file into fqdfu */
	if (opt->initialization_UMI) {
		if((err = read_initialization_file(opt->initialization_UMI,
							 &fqdfu, opt->info)))
			goto CLEAR_AND_EXIT;
		opt->K_UMI = fqdfu->n_reads;
	}

	/* create model
	 * [TODO] reasonable defaults, but uses data for binned quality models
	 */
	if ((err = make_model(&mod, dat, opt)))
		goto CLEAR_AND_EXIT;

	/* make initializer */
	if ((err = make_initializer(&ini, dat, opt,fqdf,fqdfu)))
		goto CLEAR_AND_EXIT;

	/* create run_info object */
	if ((err = make_run_info(&ri, dat, opt)))
		goto CLEAR_AND_EXIT;

	/* merely report abundances of the options::most_abundant unique
	 * sequences and exit
	 */
	if (opt->most_abundant) {
		for (unsigned int i = 0; i < opt->most_abundant; ++i)
			fprintf(stdout, " %d", ini->uniq_seq_count[i]);
		fprintf(stdout, "\n");
                goto CLEAR_AND_EXIT;
	}

	/* error estimation */
	if (opt->error_estimation) {
		if ((err = error_profile_generator(opt, dat, mod, ini, ri)))
			goto CLEAR_AND_EXIT;
		goto CLEAR_AND_EXIT;
	}

	/* main algorithm */
	if (opt->UMI_length) {
		if ((err = EM_algorithm(opt, dat, mod, ini, ri)))
			return err;

	} else if ((!opt->initialization_file) && opt->run_amplici) {

		if (opt->partition_file) {
			if ((err = ampliCI_wpartition(opt, dat, mod, ini, ri)))   // currently too slow ....
				return err;
		} else {
			if ((err = ampliCI(opt, dat, mod, ini, ri)))
				goto CLEAR_AND_EXIT;
		}

	/* reads assignment with user-provided haplotypes */
	} else if (opt->initialization_file) {
		if ((err = reads_assignment(opt, dat, mod, ini, ri)))
			goto CLEAR_AND_EXIT;

	}

CLEAR_AND_EXIT:

	if (dat)
		free_data(dat);
	if (ini)
		free_initializer(ini, opt);
	if (mod)
		free_model(mod);
	if (ri)
		free_run_info(ri);
	if (opt)
		free_options(opt);
	if(fqo)
		free_fastq_options(fqo);
	if(fqdf)
		free_fastq(fqdf);

	return (EXIT_FAILURE);	/* return err? */
} /* main */
