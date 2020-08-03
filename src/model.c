/**
 * @file model.c
 * @author Karin S. Dorman
 * @author Xiyu Peng
 *
 * Manipulate model object.
 *
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include "model.h"
#include "data.h"
#include "options.h"

/**
 * Create model object.
 *
 * @param mod	model object to create
 * @param dat	pointer to data object
 * @param opt	pointer to options object
 * @return	error status
 */
int make_model(model **mod, data *dat, options *opt)
{
	model *rm;
	*mod = malloc(sizeof **mod);

	if (*mod == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model");

	rm = *mod;
	/* K */
	rm->K = opt->K;

	/* pi */
	rm->pi = malloc(rm->K * sizeof *rm->pi);

	if (!rm->pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model.pi");

	rm->n_quality = dat->n_quality;


	/* haplotypes */
	/* [XY] Currently use ini->seeds to store haplotypes */
	rm -> haplotypes = NULL;  
	//rm->haplotypes = malloc(dat->max_read_length * opt->K
	//	* sizeof *rm->haplotypes);

	//if (!rm->haplotypes)
	//	return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
	//		"model::haplotypes");

	/* eik */
	rm->eik = malloc(dat->sample_size * rm->K * sizeof *rm->eik);

	if (rm->eik == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "model::eik");

	rm->error_profile = NULL;
	rm->err_encoding = opt->err_encoding;  

	/* compute Pr(#{indel} <= dat->max_read_length) */
	rm->adj_trunpois = ppois(dat->max_read_length, dat->max_read_length * opt->indel_error, 1, 1); //log version

	//debug_msg(DEBUG_I, DEBUG_I, "adj: %8.2e\n", rm->adj_trunpois);

 	if (opt->use_error_profile && opt->error_profile_name) {

		int fxn_debug = ABSOLUTE_SILENCE;
		rm->error_profile = calloc(NUM_NUCLEOTIDES * NUM_NUCLEOTIDES
				* rm->n_quality, sizeof *rm->error_profile);

		if (!rm->error_profile)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							 "error profile");
		
		unsigned int sq_num_nuc = NUM_NUCLEOTIDES * NUM_NUCLEOTIDES;
		double rate;

		/* may already be set at the beginning */
		//rm->err_encoding = STD_ENCODING; // usually the input error profile shows the order of A, C, G, T
		
		FILE *file = fopen(opt->error_profile_name, "rb");
		if (!file)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->error_profile_name);

		for (unsigned int r = 0; r < sq_num_nuc; r++){
			for (unsigned int q = MIN_ASCII_QUALITY_SCORE;
				 q <= MAX_ASCII_QUALITY_SCORE; q++){
				fscanf(file, "%lf,", &rate);
				if (q >= dat->fdata->min_quality && q <= dat->fdata->max_quality)
					rm->error_profile[r * rm->n_quality +
									  q - dat->fdata->min_quality] = log(rate / 1000); // log(error_rate)
				debug_msg(DEBUG_I, fxn_debug,
						  "r: %i,q:%i,rate: %f\n", r, q, rate);
			}
			debug_msg(DEBUG_I, fxn_debug, "\n");
		}
		debug_msg(DEBUG_I, fxn_debug, "finish read error profile \n");

		fclose(file);
		
	}

	/* count parameters */
	rm->n_param = opt->K * dat->max_read_length	/* haplotypes */
		+ opt->K - 1				/* pi */
		+ (opt->use_error_profile?dat->n_quality*12:0);
	
	rm->ll = -INFINITY;
	rm->best_ll = -INFINITY;	// best log likelihood
	rm->JC_ll = -INFINITY;

	rm->aic = INFINITY;
	rm->bic = INFINITY;

	rm->distance = NULL;
	rm->est_ancestor = NULL;
	rm->JC_ll_K = NULL;

	if (opt->JC69_model) {
		rm->distance = malloc(rm->K * sizeof *rm->distance);
		rm->est_ancestor = malloc(dat->max_read_length
					* sizeof * rm->est_ancestor);
		rm->JC_ll_K = malloc(rm->K * sizeof *rm->JC_ll_K);

		if (!rm->distance || !rm->est_ancestor || !rm->JC_ll_K)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"model::JC69");
	}

	return NO_ERROR;
} /* make_model */

/**
 * Reallocate model struct for a different K or different sample with different size.
 *
 * @param mod	pointer to model object
 * @param dat	pointer to data object
 * @param opt	pointer to options objet
 * @return	error status
 */
int realloc_model(model *mod, data *dat, options *opt)
{

	mod->K = opt->K;

	/* pi  (store in log(pi)) */
	double *pi = realloc(mod->pi, mod->K * sizeof *mod->pi);
	if (!pi)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.model.pi");
	mod->pi = pi;


	/* haplotypes */
	if (!dat->max_read_position || !mod->K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"realloc.model.haplotypes");
	/* 
	unsigned char *haplotypes = realloc(mod->haplotypes,
		dat->max_read_length * mod->K * sizeof *mod->haplotypes);

	if (!haplotypes)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
						"realloc.model.haplotypes");
	mod->haplotypes = haplotypes;
	*/

	/* e_ik */
	if (!dat->sample_size || !mod->K)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "realloc.model.eik");
	double *eik = realloc(mod->eik,
		dat->sample_size * mod->K * sizeof *mod->eik);
	if (!eik)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "realloc.model.eik");
	mod->eik = eik;

	mod->best_ll = -INFINITY;	/* best log likelihood */


	//mod->n_quality = dat->n_quality;
	/* count parameters */
	mod->n_param = opt->K * dat->max_read_length	/* haplotypes */
		+ opt->K - 1				/* pi */
		+ (opt->use_error_profile?dat->n_quality*12:0);  
		//+ 12					/* gamma */
		//+ (opt->background_model ? NUM_NUCLEOTIDES - 1 : 0);	/* bg_pi */


	mod->aic = INFINITY;
	mod->bic = INFINITY;

	if(opt->JC69_model){
		mod->distance = realloc(mod->distance, mod->K * sizeof *mod->distance);
		mod->JC_ll_K = realloc(mod->JC_ll_K, mod->K * sizeof *mod->JC_ll_K);
	}

	mmessage(INFO_MSG, NO_ERROR, "Number of parameters: %u\n", mod->n_param);

	return NO_ERROR;
}/* realloc_model */


/**
 * Extract error probability from error profile. (log version)
 *
 * [TODO] Consider making this an inline function.
 * @param error_profile	the error profile
 * @param n_quality	the number of possible quality scores
 * @param hap_nuc	true nucleotide
 * @param obser_nuc	observed nucleotide (with error, possibility)
 * @param qual		observed quality score
 * @return		error probability
 */
double translate_error_STD_to_XY(double *error_profile, unsigned char n_quality,
	unsigned char hap_nuc, unsigned char obser_nuc, unsigned char qual){

	double lp = 0.;

     //A (0),C(1),G(3),T(2)
	if (hap_nuc == XY_A) {
		if (obser_nuc == XY_A)
			lp = error_profile[qual]; //A->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality+qual];  //A->C
		else if (obser_nuc == XY_T)   //A->T
			lp = error_profile[n_quality*3+qual];
		else  // A->G
			lp = error_profile[n_quality*2+qual];
	} else if (hap_nuc == XY_C) {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality*4+qual]; //C->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality*5+qual];  //C->C
		else if (obser_nuc == XY_T)   //C->T
			lp = error_profile[n_quality*7+qual];
		else  // C->G
			lp = error_profile[n_quality*6+qual];
	} else if (hap_nuc == XY_G) {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality*8+qual]; //G->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality*9+qual];  //G->C
		else if (obser_nuc == XY_T)   //G->T
			lp = error_profile[n_quality*11+qual];
		else  // G->G
			lp = error_profile[n_quality*10+qual];
	} else {
		if (obser_nuc == XY_A)
			lp = error_profile[n_quality*12+qual]; //T->A
		else if (obser_nuc == XY_C)
			lp = error_profile[n_quality*13+qual];  //T->C
		else if (obser_nuc == XY_T)   //T->T
			lp = error_profile[n_quality*15+qual];
		else  // T->G
			lp = error_profile[n_quality*14+qual];
	}

	return lp;
}/* translate_error_STD_to_XY */

/**
 * Delete model object.
 *
 * @param mod	pointer to model object to delete
 */
void free_model(model *mod)
{
	if (!mod)
		return;
	if (mod->pi) free(mod->pi);
	if (mod->eik) free(mod->eik);
	if (mod->distance) free(mod->distance);
	if (mod->JC_ll_K) free(mod->JC_ll_K);
	if (mod->haplotypes) free(mod->haplotypes);
	if (mod->est_ancestor) free(mod->est_ancestor);
	if (mod->error_profile) free(mod->error_profile);
	free(mod);
	mod = NULL;
} /* free_model */
