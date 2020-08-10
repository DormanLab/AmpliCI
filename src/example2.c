/**
 * @file example2.c
 * @author Xiyu Peng
 *
 * An example for calling amplici_wfile function
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include "libamplici.h"
#include <stdio.h>
#include <stdarg.h>


int main(int argc, const char **argv)
{
    UNUSED(argc);
    UNUSED(argv);

    int err = 0;

    /* read two input files */
    char *fastq_file = "../test/SRR2990088_1_noN_3000_subset1_new.fastq";
    char *error_profile_name = NULL;
    char *output_file = "test.out";

    /* initialize output */
    unsigned int K = 0;
    unsigned int *cluster_id = NULL;
    unsigned int *cluster_size = NULL;
    unsigned char *seeds = NULL;
    unsigned int *seeds_length = NULL;
    size_t sample_size = 0;

    if ((err = amplici_wfile(fastq_file, error_profile_name, &seeds, &seeds_length, &cluster_id,
                       &cluster_size, &K, &sample_size)))
        goto CLEAR_AND_EXIT;

    /* print and check */
    FILE *fp = NULL;
    fp = fopen(output_file, "w");
    if (!fp)
        return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
                        output_file);

    fprintf(fp, "K: %i\n", K);

    fprintf(fp, "assignments: ");
    fprint_assignment(fp, cluster_id, sample_size,
                      K, 2, 1);
    fprintf(fp, "cluster sizes: ");

    fprint_uints(fp, cluster_size, K, 3, 1);
    
    fprint_fasta(fp, seeds, K, seeds_length[0],seeds_length, "H");
    
    fclose(fp);

CLEAR_AND_EXIT:
    if (cluster_id)
        free(cluster_id);
    if (cluster_size)
        free(cluster_size);
    if (seeds)
        free(seeds);
    if (seeds_length)
        free(seeds_length);

    return (EXIT_FAILURE);
} /* main */
