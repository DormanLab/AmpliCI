/**
 * @file partition.h
 * @author Xiyu Peng
 *
 * Header Files for partition.c
 */

#ifndef __H_PARTITION__
#define __H_PARTITION__

#include "amplici.h"

int read_partition_file(char const * const filename, unsigned int * cluster_id, unsigned int sample_size);
int ampliCI_wpartition(options * opt, data * dat, model *mod, initializer *ini, run_info *ri);

#endif

