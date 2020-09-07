/**
 * @file example.c
 * @author Xiyu Peng
 *
 * An example for calling amplici_core function
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include "libamplici.h"
#include "fastq.h"
#include "error.h"
#include "amplici.h"
#include "io.h"

int main()
{   

    int err = NO_ERROR;

    /* read two input files */
    char *fastq_file = "./test/SRR2990088_1_noN_3000_subset1_new.fastq";
    char *error_profile_name =   NULL;  
	char *output_file = "test.out";
    double low_bound = 2.0;

    /* initialize output */
    unsigned int K = 0;
    unsigned int *cluster_id = NULL;
    unsigned int *cluster_size = NULL;
    unsigned char *seeds = NULL;
	unsigned int *seeds_length = NULL;
    double *ll = NULL;
    double *abun = NULL;

    /* initialize input */
    fastq_options *fqo = NULL;	/* fastq file options */
    fastq_data *fdata = NULL;
    unsigned char **dmat = NULL;
    unsigned char **qmat = NULL;

    if ((err = make_fastq_options(&fqo)))
		goto CLEAR_AND_EXIT;


	/* encode nucleotides in 2-bits: error raised if ambiguous bases */
	fqo->read_encoding = XY_ENCODING;



	/* read sequence data */
	if (!fastq_file || (err = read_fastq(fastq_file,
						&fdata, fqo)))
		goto CLEAR_AND_EXIT;

    size_t sample_size = fdata->n_reads;
    unsigned int max_read_length = fdata->n_max_length;

    dmat = malloc(sample_size * sizeof *dmat);

	if (!dmat)
		goto CLEAR_AND_EXIT;

	unsigned char *rptr = fdata->reads;
	for (size_t i = 0; i < sample_size; ++i) {
		dmat[i] = rptr;
		rptr += max_read_length;
	}

	/* allocate short-cut pointers to quality sequences  */
	qmat = malloc(sample_size * sizeof *qmat);

	if (!qmat)
		goto CLEAR_AND_EXIT;

	unsigned char *qptr = fdata->quals;
	for (size_t i = 0; i < sample_size; i++) {
		qmat[i] = qptr;
		qptr += max_read_length;
	}


    /* amplicon clustering */
    if((amplici_core(dmat, qmat, sample_size, low_bound, max_read_length, 
                error_profile_name, fdata->max_quality, fdata->min_quality, 
                &seeds, &seeds_length, &cluster_id, 
                &cluster_size, &K, &ll, &abun)))
        goto CLEAR_AND_EXIT;

   
    /* print and check */
	FILE *fp = NULL;
    fp = fopen(output_file, "w");
    if (!fp)
        goto CLEAR_AND_EXIT;

    fprintf(fp, "K: %i\n", K);

    fprintf(fp, "assignments: ");
    fprint_assignment(fp, cluster_id, sample_size,
                      K, 2, 1);
    fprintf(fp, "cluster sizes: ");

    fprint_uints(fp, cluster_size, K, 3, 1);
    
    fprint_fasta(fp, seeds, K, max_read_length, seeds_length, "H");
    
    fclose(fp);


CLEAR_AND_EXIT:
    if(fdata) free_fastq(fdata);
    if(fqo) free_fastq_options(fqo);
    if (dmat) free(dmat);
	if (qmat) free(qmat);
	if (cluster_id) free(cluster_id);
	if (cluster_size) free(cluster_size);
	if (seeds)  free(seeds);		
	if (seeds_length) free(seeds_length);
    if (ll) free(ll);
    if (abun) free(abun);


    return(EXIT_FAILURE);
}/* main */