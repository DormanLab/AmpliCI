/**
 * @file amplici_umi.h
 * @author Xiyu Peng
 *
 * Cluster amplicon sequences with Unique Molecular Identifier
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include "options.h"
#include "initialize.h"
#include "amplici.h"

int ini_eta_gamma(options *opt, initializer *ini, double *gamma, double *eta, size_t sample_size);
int EM_algorithm(options *opt, data *dat, model *mod, initializer *ini, run_info *ri);