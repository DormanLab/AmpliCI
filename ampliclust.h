/**
 * @file ampliclust.h
 * @author Karin S. Dorman
 *
 * Header file for ampliclust.
 */

#ifndef __H_AMPLICLUST__
#define __H_AMPLICLUST__

#include <stdlib.h>
#ifdef USE_CURSES
#include <curses.h>
#endif

#include "model.h"
#include "data.h"
#include "options.h"
#include "initialize.h"


/**
 * Where to cache results.  Is this the current solution, the best local
 * solution, the optimal global solution, or the truth?
 */
enum {
	DEFAULT,	/*<! cluster_id, criterion, etc. */
	BEST,		/*<! best_cluster_id, best_criterion, etc. */
	OPTIMAL,	/*<! optimal_cluster_id, etc. */
	TRUE_VALUES,	/*<! true values */
	INITIAL_VALUES	/*<! initial values */
};

typedef struct _run_info run_info;

/**
 * Store information about the best solution across multiple runs (varying 
 * initialization or K or whatever).
 */
struct _run_info {
	
	unsigned int *optimal_cluster_id;	/*<! optimal hard clustering */
	unsigned int *optimal_cluster_size;	/*<! optimal cluster sizes */
	double *optimal_cluster_ll; /*<! maximal read log likelihood */

};

int assign_clusters(double *eik, unsigned int K, size_t sample_size, unsigned int *cs, unsigned int *ci,int by_K); 
void fprint_fasta(FILE *fp, data_t *data, size_t n, size_t p, char const * const prefix);  
void fprint_alignment2(FILE *fp, data_t **data, size_t n, size_t p);
#ifdef USE_CURSES
void wprint_fasta(WINDOW *wp, data_t *data, size_t n, size_t p, char const * const prefix);
void wprint_alignment2(WINDOW *wp, data_t **data, size_t n, size_t p);
#endif

int make_run_info(run_info **ri, data *dat,options *opt);
int realloc_run_info(run_info *ri,size_t sample_size,unsigned int K,int change_sample, int change_K);
void free_run_info(run_info *ri);

#endif
