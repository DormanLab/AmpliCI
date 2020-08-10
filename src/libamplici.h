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
#include "fastq.h"
#include "amplici.h"
#include "io.h"

int amplici_core(data_t **dmat, data_t **qmat, size_t sample_size, unsigned int rlen, 
                double *error_profile, unsigned int n_quality, data_t ***seeds, unsigned int **seeds_length, unsigned int **cluster_id, 
                unsigned int** cluster_size, unsigned int *K);

int amplici_wfile(char *fastq_file, char *error_profile_name, unsigned char **seeds, unsigned int **seeds_length, 
               unsigned int **cluster_id, unsigned int **cluster_size, unsigned int *K, size_t *sample_size);

