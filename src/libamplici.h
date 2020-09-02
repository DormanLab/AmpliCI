/**
 * @file libamplici.h
 * @author Xiyu Peng
 *
 * Header file for calling amplici_core function
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include "constants.h"

int amplici_wfile(char *fastq_file, char *error_profile_name, double low_bound, unsigned char **seeds, unsigned int **seeds_length,
                  unsigned int **cluster_id, unsigned int **cluster_size, unsigned int *K, size_t *sample_size,
                  unsigned int *max_read_length, double **abun, double **ll);

int amplici_core(data_t **dmat, data_t **qmat, size_t sample_size, unsigned int max_read_length,
                 char *error_profile_name, unsigned int n_quality, unsigned char min_quality, unsigned char **seeds, unsigned int **seeds_length,
                 unsigned int **cluster_id, unsigned int **cluster_size, unsigned int *K);
