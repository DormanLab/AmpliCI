/**
 * @file amplici_umi.c
 * @author Xiyu Peng
 *
 * Cluster amplicon sequences with Unique Molecular Identifier
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 * */


#include "options.h"
#include "initialize.h"
#include "amplici.h"
#include "align.h"
#include "model.h"
#include "io.h"
#include "amplici_umi.h"
#include "error.h"

#define MATHLIB_STANDALONE 1
#include <Rmath.h>


int trans_hap_and_umi(options *opt, data *dat, initializer *ini, model *mod);
double E_step(model *mod, size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI, int *err);
int M_step(options* opt, model *mod,size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI);
int reads_assign_sparse(model *mod, run_info *ri, size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI);
int reads_assign_optimal(model *mod, run_info *ri, size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI);
double MPLE_gamma_s(double *x_s, unsigned int K, int *err, unsigned int s,double rho, double omega);
int trans_expect_UMIs(options *opt, data *dat, data_t *seeds_UMI, double *error_profile,
						double *trans_prob, int ends_free);


/* main model for UMI */
int EM_algorithm(options *opt, data *dat, model *mod, initializer *ini, run_info *ri)
{
	int err = NO_ERROR;
	int output_hap = 1;
	double rdelta;
	unsigned int EM_iter;
	double epsilon = opt->epsilon;

	/* calculate the transition prob */
	if ((err = trans_hap_and_umi(opt, dat, ini, mod)))
		return err;


	/* initialize gamma and eta */
	mod->penalty_ll = ini_eta_gamma(opt, ini, mod->gamma, mod->eta, dat->sample_size);
	if (isnan(mod->penalty_ll))
		fprintf(stderr, "penalty_ll : %15.7f \n", mod->penalty_ll);

	/* when ll relative diff > 1e-6 */

	/* simple E step, output E_id and E_value and ll given pgamma and peta */

	EM_iter = 0;

	do {

		mmessage(INFO_MSG, NO_ERROR, "The E step of the %5u iteration...\n", EM_iter);

		// print gamma
		/*
		fprintf(stderr, "\nEst gamma: ");
		for (unsigned int k = 0; k < opt->K; ++k) {
			fprintf(stderr, "%15.3f\t", mod->gamma[326*opt->K+k]);
		}
		*/

		log_vector(mod->gamma, opt->K*opt->K_UMI);
		log_vector(mod->eta, opt->K_UMI);

		mod->ll_UMI = E_step(mod, dat->sample_size, opt->topN,
					opt->K, opt->K_UMI, &err);
		mod->ll_UMI -= mod->penalty_ll;
		if(err)
			return err;


		rdelta = (mod->pll_UMI - mod->ll_UMI)/mod->pll_UMI;
		mmessage(INFO_MSG, NO_ERROR, "%5u: %15.3f (%15.3f - %15.3f, %5.3e)\n", EM_iter,
					mod->ll_UMI,mod->ll_UMI+mod->penalty_ll, mod->penalty_ll, rdelta);


		mmessage(INFO_MSG, NO_ERROR, "The M step of the %5u iteration...\n", EM_iter);

		/* update gamma and eta */
		if ((err = M_step(opt,mod,dat->sample_size, opt->topN,
					opt->K, opt->K_UMI)))
			return err;


		mod->pll_UMI = mod->ll_UMI;
		EM_iter ++;

		if (rdelta < epsilon && rdelta > 0)
			break;


		if (EM_iter > opt->n_iter_amplici) {
			mmessage(WARNING_MSG, INTERNAL_ERROR,
				"exceed max interations,%u\n", EM_iter);
			break;
		}

	} while(1);

	mmessage(INFO_MSG, NO_ERROR, "Final reads assignment.....\n");

	/* reads assignment */

	//if((err = reads_assign_sparse(mod, ri, dat->sample_size, opt->topN,
	//				opt->K, opt->K_UMI)))
	//	return err;


	/* assignment with the combination with the joint maximum log likelihood */
	reads_assign_optimal(mod, ri, dat->sample_size, opt->topN,
				   opt->K, opt->K_UMI);


	/* filtering on gamma based on the threshold on eta [Try remove false positive UMIs ] */
	double thres_eta = opt->threshold_UMI / dat->sample_size;
	for (unsigned int b = 0; b < opt->K_UMI; ++b) {
		if (mod->eta[b] < thres_eta) {  // filter UMIs with abundance < the threshold
			for (unsigned int k = 0; k < opt->K; ++k)
				mod->gamma[b*opt->K + k] = 0.;
		   // fprintf(stderr, "remove the %d th UMI\n",b);
		}
	}


	if (opt->umicollision) { // allow collision
		for (unsigned int k = 0; k < opt->K; ++k) {
			for (unsigned int b = 0; b < opt->K_UMI; ++b) {
				// post hoc separate UMI+ haplotypes if there is an UMI collision
				if (mod->gamma[b*opt->K + k] > 0 && mod->eta[b] > thres_eta)
					mod->gamma[b*opt->K + k] = 1.0;
			}
		}
	} else { // Or just assign UMIs to the haplotype with the highest probability
		for (unsigned int b = 0; b < opt->K_UMI; ++b) {
			double max = 0.;
			unsigned int max_idx = 0;

			if (mod->eta[b] < thres_eta)
				continue;

			for(unsigned int k = 0; k < opt->K; ++k) {
				if (mod->gamma[b*opt->K + k] > max) {
					max = mod->gamma[b*opt->K + k];
					max_idx = k;
				}
			}
			for (unsigned int k = 0; k < opt->K; ++k) {
				if (k == max_idx)
					mod->gamma[b*opt->K + k] = 1.0;
				else
					mod->gamma[b*opt->K + k] = 0.;
			 }
		}
	}

	/* Deduplicate abundance for each haplotypes */
	for (unsigned int k = 0; k < opt->K; ++k) {
		ini->H_abun[k] = 0.;
		for (unsigned int b = 0; b < opt->K_UMI; ++b)
			ini->H_abun[k] += mod->gamma[b*opt->K + k];
	}

	char *outfile_hap = NULL;
	char *outfile = NULL;

	if (opt->outfile_base || opt->outfile_info) {
		/* print result and check */
		FILE *fp = NULL;

		if (!opt->outfile_info) {
			outfile = malloc((strlen(opt->outfile_base) + 5) * sizeof(char));
			if (!outfile)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"output file");
			strcpy(outfile, opt->outfile_base);
			strcat(outfile, ".out");
			opt->outfile_info = outfile;
		}

		fp = fopen(opt->outfile_info, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->outfile_info);

		fprintf(fp, "log likelihood: %f\n", mod->ll_UMI);

		fprintf(fp, "K: %i\n", opt->K);

		fprintf(fp, "read assignments: ");
		fprint_assignment(fp, ri->optimal_cluster_id, dat->sample_size,
						  opt->K, 2, 1);

		fprintf(fp, "cluster sizes: ");

		fprint_uints(fp, ri->optimal_cluster_size, opt->K, 3, 1);

		fprintf(fp, "UMI K: %i\n", opt->K_UMI);

		fprintf(fp, "UMI assignments: ");
		fprint_assignment(fp, ri->UMI_cluster_id, dat->sample_size,
						  opt->K_UMI, 2, 1);

		fprintf(fp, "UMI cluster sizes: ");

		fprint_uints(fp, ri->UMI_cluster_size, opt->K_UMI, 3, 1);

		fprintf(fp, "reads ll: ");
		fprint_doubles(fp, ri->optimal_cluster_ll, dat->sample_size, 3, 1);

		fprintf(fp, "Eta: ");
		fprint_doubles(fp, mod->eta, opt->K_UMI, 6, 1);

		fprintf(fp, "Gamma: ");
		fprint_vectorized_matrix(fp, mod->gamma, opt->K_UMI, opt->K, 1);

		/* output haplotypes and related UMIs */
		unsigned int cnt = 0;
		unsigned int l = opt->UMI_length;
		for (unsigned int k = 0; k < opt->K; ++k) {
			if (ini->H_abun[k] > opt->threshold_hap) {
				for (unsigned int j = 0; j < ini->seed_lengths[k]; ++j)
					fprintf(fp, "%c", xy_to_char[(int)ini->seeds[k][j]]);
				fprintf(fp, "\n");
				for (unsigned int b = 0; b < opt->K_UMI; ++b) {
					if(mod->gamma[b*opt->K + k] > 0.) {
						for (unsigned int j = 0; j < opt->UMI_length; ++j)
						 fprintf(fp, "%c", xy_to_char[(int)ini->seeds_UMI[b * l +j]]);
						fprintf(fp, "\t");
					}
				 }
				 fprintf(fp, "\n");
				cnt++;
			}
		}

		fclose(fp);

		mmessage(INFO_MSG, NO_ERROR, "Output the final result file: "
				 "%s \n", opt->outfile_info);
	}
	if(outfile)
		free(outfile);

	if (output_hap) {

		FILE *fp2 = NULL;

		if (!opt->outfile_fasta) {
			if (!opt->outfile_base)
				return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"invalid output filenames");

			outfile_hap = malloc((strlen(opt->outfile_base)
							+ 5) * sizeof(char));
			if (!outfile_hap)
				return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
								"output file");

			strcpy(outfile_hap, opt->outfile_base);
			strcat(outfile_hap, ".fa");
			opt->outfile_fasta = outfile_hap;
		}

		fp2 = fopen(opt->outfile_fasta, "w");
		if (!fp2)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
							opt->outfile_fasta);

		unsigned int cnt = 0;
		//char const * const prefix = "H";
		for (unsigned int k = 0; k < opt->K; ++k) {
			if (ini->H_abun[k] > opt->threshold_hap) {
				fprintf(fp2, ">%s%d;", "H", cnt);
				fprintf(fp2,"Deduplicated Abundance=%.3f;",ini->H_abun[k]);
				fprintf(fp2, "\n");
				for (unsigned int j = 0; j < ini->seed_lengths[k]; ++j)
					fprintf(fp2, "%c", xy_to_char[(int)ini->seeds[k][j]]);
				fprintf(fp2, "\n");
				cnt++;
			}
		}

		fclose(fp2);

		mmessage(INFO_MSG, NO_ERROR, "Output the final haplotype fasta "
					"file: %s \n", opt->outfile_fasta);
	}

	if (outfile_hap)
		free(outfile_hap);

	return err;
}/* EM_algorithm */


/* set initial value for gamma and eta */
double ini_eta_gamma(options *opt, initializer *ini, double *gamma, double *eta, size_t sample_size)
{

	//int err = NO_ERROR;
	unsigned int K = opt->K;
	int fxn_debug = ABSOLUTE_SILENCE;
	double llp_sum;
	double rho = opt->rho;
	double omega = opt->omega;


	for (size_t i = 0; i < sample_size; i++) {
		if (ini->reads_umi_id[i] > -1) {
			eta[ini->reads_umi_id[i]] += 1.0;
			if (ini->reads_hap_id[i] > -1)
				gamma[ini->reads_umi_id[i] * K + ini->reads_hap_id[i]] += 1.0;
		}
	}

	double eta_sum = 0.;
	unsigned int idx;

	llp_sum = 0.;
	for (unsigned int b = 0; b < opt->K_UMI; ++b) {
		eta_sum += eta[b];

		debug_msg(DEBUG_II, fxn_debug, "%15.3f,",eta[b]);
		debug_msg(DEBUG_II, fxn_debug, "\n");

		double gamma_sumh = 0.;
		for (unsigned int k = 0; k < K; ++k)
			gamma_sumh += gamma[b * K + k];

		for (unsigned int k = 0; k < K; ++k) {
			idx = b * K + k;
			debug_msg(DEBUG_II, fxn_debug, "%15.3f / %15.3f\n,",gamma[idx],gamma_sumh);
			gamma[idx] /=  (gamma_sumh+ 1e-40); // avoid /0
			debug_msg(DEBUG_III, fxn_debug, "%15.3f,",gamma[idx]);
			double llp = rho * log1p(gamma[idx]/omega);
			llp_sum += llp;
			if (isnan(llp))
				fprintf(stderr, "(%d, %d) : gamma = %15.7f, ll = %15.7f \n", b, k, gamma[idx], llp);
		}
		debug_msg(DEBUG_III, fxn_debug, "\n");
	}

	for (unsigned int b = 0; b < opt->K_UMI; ++b)
		eta[b] /= eta_sum;

	if (opt->trans_penalty == MLE)
		return 0;

	return llp_sum;
}/* ini_eta_gamma */

/* calculate transition prob between data and input UMI and Hap */
int trans_hap_and_umi(options *opt, data *dat, initializer *ini, model *mod)
{

	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;

	double *error_profile = NULL;
	if (opt->use_error_profile && mod->error_profile) {
		error_profile = mod->error_profile;
		debug_msg(DEBUG_II, fxn_debug, "Use error profile. \n");
	}

	mmessage(INFO_MSG, NO_ERROR, "Calculating the trans prob of reads......\n");


	/* For acceleration purpose. Do not allow large gap in amplification and sequencing */
	opt->band = 9;	// still need check
	opt->ends_free = 1;

	/* transition probability for reads */
	if ((err = trans_expectation(opt, dat, ini, error_profile,
			 mod->adj_trunpois, mod->eik, opt->ends_free)))
		return err;


	/* transition probability for UMIs */
	/* implement a naive method, waiting for an efficient one */

	mmessage(INFO_MSG, NO_ERROR, "Calculating the trans prob of UMIs......\n");

	 /* [TODO] Set a new model of UMI alignment, not use nw alignment */

	opt->band = 0;   // for UMIs // need further investigation
	// opt->gap_p = -20;
	// opt->ends_free = 1;

	/* transition probability for UMIs */
	if ((err = trans_expect_UMIs(opt, dat, ini->seeds_UMI, error_profile, mod->eik_umi,opt->ends_free)))
		return err;
	/* calculate the transition probability without alignment */

	/*
	normalize(dat->sample_size, opt->K, mod->eik);
	normalize(dat->sample_size, opt->K_UMI, mod->eik_umi);
	*/

	return err;
}/* trans_hap_and_umi */

/* update E2_sparse */
double E_step(model *mod, size_t sample_size, unsigned int topN,
			unsigned int K, unsigned int K_UMI, int *err)
{

	double ll = 0.;
	unsigned int idx;
	unsigned int s;
	double max, cur_value;
	double maximum;
	double sum;
	int fxn_debug = ABSOLUTE_SILENCE;
	double read_ll;
	// double log_epsilon = - 40;

	/* [TODO] temp is sparse, think about a more efficient way to do */
	double *temp = NULL;
	temp = calloc(K* K_UMI, sizeof *temp);
	if (!temp) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "temp");
		return -INFINITY;
	}

	double *gamma = mod->gamma;
	double *eta = mod->eta;


	/* Find top N log joint probabilities for each read */
	for (unsigned int i = 0; i < sample_size; ++i) {
		idx = i * topN;
		s = 0;

		/* Calculation and find the maximum value */
		max = -INFINITY;
		for (unsigned int k = 0; k < K; ++k) {
			for (unsigned int b = 0; b < K_UMI; ++b){
				if (gamma[b*K+k] > -INFINITY && eta[b] > -INFINITY) {
					cur_value = gamma[b*K + k] + eta[b]
						+ mod->eik_umi[b*sample_size +i] + mod->eik[k*sample_size + i];
				} else {
					debug_msg(DEBUG_I, fxn_debug, "gamma == 0\n");
					cur_value = -INFINITY;
				}
				temp[b*K+k] = cur_value;
				if (cur_value > max) {
					max = cur_value;
					mod->E2_sparse_value[idx + s] = cur_value;
					mod->E2_sparse_hap_id[idx + s] = k;
					mod->E2_sparse_umi_id[idx + s] = b;
				}
			}
		}
		if (i < 10) {
			debug_msg(DEBUG_I, fxn_debug, "maximum value for the %5d reads: %15.3f\n",i, max);
			debug_msg(DEBUG_I, fxn_debug, "gamma : %15.3f\n",gamma[mod->E2_sparse_umi_id[idx + s]*K+mod->E2_sparse_hap_id[idx + s]]);
			debug_msg(DEBUG_I, fxn_debug, "eta : %15.3f\n",eta[mod->E2_sparse_umi_id[idx + s]]);
			debug_msg(DEBUG_I, fxn_debug, "EUMI : %15.3f\n",mod->eik_umi[mod->E2_sparse_umi_id[idx + s]*sample_size +i]);
			debug_msg(DEBUG_I, fxn_debug, "E : %15.3f\n",mod->eik[mod->E2_sparse_hap_id[idx + s]*sample_size +i]);
		}

		temp[mod->E2_sparse_umi_id[idx + s]*K+mod->E2_sparse_hap_id[idx + s]] = -INFINITY;
		s += 1;
		maximum = max;  // maximum value among top N

		while (s < topN) {  //Need a more efficient algorithm here
			max = -INFINITY;
			for (unsigned int k = 0; k < K; ++k) {
				for (unsigned int b = 0; b < K_UMI; ++b) {
					cur_value = temp[b * K + k];
					if (cur_value > max) {
						max = cur_value;
						mod->E2_sparse_value[idx + s] = cur_value;
						mod->E2_sparse_hap_id[idx + s] = k;
						mod->E2_sparse_umi_id[idx + s] = b;
					}
				}
			}
			temp[mod->E2_sparse_umi_id[idx + s]*K+mod->E2_sparse_hap_id[idx + s]] = -INFINITY;
			s += 1;

			debug_msg(DEBUG_II, fxn_debug, "maximum value for the %5d reads: %15.3f\n",i, max);

		}
		sum = 0.;
		for (unsigned int t = 0; t < topN; ++t) {
			mod->E2_sparse_value[idx + t] = exp(mod->E2_sparse_value[idx + t] - maximum);
			sum += mod->E2_sparse_value[idx + t];
		}

		read_ll = log(sum) +maximum;
		ll += read_ll;
		if (i < 10)
			debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, read_ll);

		for (unsigned int t = 0; t < topN; ++t)
			mod->E2_sparse_value[idx + t] /= sum;

	}
	debug_msg(DEBUG_III, fxn_debug, "\n");

	if (temp)
		free(temp);

	return ll;
}/* E_step */



/* estimate new eta and new gamma */
int M_step(options *opt, model *mod,size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI)
{
	int err = NO_ERROR;
	//int fxn_debug = DEBUG_I;
	unsigned int idx;
	double epsilon = 1e-40;
	int algo = opt->trans_penalty;

	double *eta =  mod->eta;
	double *gamma = mod->gamma;


	for (unsigned int b = 0; b < K_UMI; ++b) {
		eta[b] = 0.;
		for (unsigned int k = 0; k < K; ++k)
			gamma[b * K + k] = 0.;
	}

	for (unsigned int i= 0; i< sample_size; ++i) {
		for (unsigned int t = 0; t < topN;++t) {
			idx = i *topN + t;
			gamma[mod->E2_sparse_umi_id[idx]* K + mod->E2_sparse_hap_id[idx]] += mod->E2_sparse_value[idx];
			//if(mod->E2_sparse_umi_id[idx] == 326 && mod->E2_sparse_hap_id[idx] == 189)
			//	fprintf(stderr, "%d: %d: %15.7f\t",i, mod->E2_sparse_hap_id[idx], mod->E2_sparse_value[idx]);
		}
	}

	for (unsigned int b = 0; b < K_UMI; ++b)
		 for (unsigned int k = 0; k < K; ++k)
			 eta[b] += gamma[b * K + k];

	//fprintf(stderr, "Est gamma: \n");
	//for(unsigned int k = 0; k < opt->K; ++k){
	//	if(mod->gamma[326*opt->K+k]>0)
	//		fprintf(stderr, "%15.3f\t", mod->gamma[326*opt->K+k]);
	//}

	//for(unsigned int k = 0; k < K; ++k)
	//   fprintf(stderr, "Est gamma: %15.3f", gamma[k]);

	// debug_msg(DEBUG_I, fxn_debug, "\n");

	/* gamma */
	mod->penalty_ll = 0.;
	for (unsigned int b = 0; b < K_UMI; b++) {

		if (algo == MPLE) {
			/* for each row of gamma */
			double lls = MPLE_gamma_s(gamma, K, &err, b,opt->rho, opt->omega);
			if(err) return err;

			/* below have beed move to the MPLE_gamma_s itself*/
			/*
			if (isnan(lls)) {
				debug_msg(DEBUG_I, fxn_debug, "ll penalty: %15.3f\n",lls);

				for(unsigned int k = 0; k < K; ++k)
					fprintf(stderr, "Est gamma: %15.3f", gamma[b * K + k]);

				fprintf(stderr, "\n");

				// compute MLE instead [XY: is it fine ? ]
				for (unsigned int k = 0; k < K; k++) {
					idx = b * K + k;
					gamma[idx] /= (eta[b] + epsilon);  // avoid /0
					double llp = opt->rho * log1p(gamma[idx]/opt->omega);
					mod->penalty_ll += llp;
				}
			*/

			//} else {
				mod->penalty_ll += lls;
			//}

		} else {
			for (unsigned int k = 0; k < K; k++) {
				idx = b * K + k;
				gamma[idx] /= (eta[b] + epsilon); // avoid /0
			}
		}
	}


	//for(unsigned int k = 0; k < K; ++k){
	//	fprintf(stderr, "Est gamma: %15.3f", gamma[k]);
	//}
	//debug_msg(DEBUG_I, fxn_debug, "\n");


	/* eta */
	for (unsigned int b = 0; b < K_UMI; b++)
		eta[b] /= sample_size;

	return err;
}/* M_step */

/* Caculate the MPLE for gamma for each s and return the panalty turn */
double MPLE_gamma_s(double *x_s, unsigned int K, int *err, unsigned int s, double rho, double omega)
{

	int fxn_debug = ABSOLUTE_SILENCE;
	double lls;
	int start = s*K;
	double max_xi = 0.;
	unsigned int arg_max_idx = 0;
	double sum_xi = 0.;

	/* set up */
	mstep_data md = {
		.omega =omega,
		.rho = rho,
		.xi = NULL,
		.idx = NULL,
		.signs = NULL,
		.n_xi = 0,
	};


	/* count the no-zeros number */
	for (unsigned int k =0 ; k < K; ++k) {
		if (x_s[start+k] > 0) {
			md.n_xi+= 1;
			sum_xi += x_s[start+k];
		}
	}

	if (md.n_xi == 0)
		return 0;

	if (md.n_xi == 1) {
		for (unsigned int k =0 ; k < K; ++k)
			if (x_s[start+k] > 0)
				x_s[start+k] = 1.0;
		return md.rho * log1p(1.0/md.omega);
	}

	md.xi = calloc(md.n_xi, sizeof (*md.xi));
	md.idx = calloc(md.n_xi, sizeof (*md.idx));
	md.signs = calloc(md.n_xi, sizeof(*md.signs));

	if (!md.xi || !md.idx || !md.signs) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err,"md\n");
		return 0;
	}

	unsigned int cnt = 0;
	for (unsigned int k = 0; k < K; ++k) {
		if (x_s[start+k] > 0) {
			 md.xi[cnt] = x_s[start+k];
			md.idx[cnt] = k;
			md.signs[cnt] = 1;
			if (md.xi[cnt]  > max_xi) {
				max_xi = md.xi[cnt];
				arg_max_idx = cnt;
			}
			cnt ++;
		}
	}

	msp_lcstr_func _plfunc = mstep_pen1_lagrange_cstr_func;
	msp_lcstr_derv _plderv = mstep_pen1_lagrange_cstr_derv;


	double l0 = 0.;
	// check the signs of the solution
	if (max_xi < md.rho) { // negative lambda

		//double lamb_lb = mstep_pen1_lambda_support(&md);
		//debug_msg(DEBUG_I, fxn_debug, "low bound for lambda: %15.3f\n",lamb_lb);
		//l0 = lamb_lb + 1e-6;

		l0 = max_xi - md.rho + 1e-6;

		// try different signs
		double l1 = l0 - mstep_pen1_lagrange_cstr_func(l0, &md) /
			mstep_pen1_lagrange_cstr_derv(l0, &md);

		md.signs[arg_max_idx] = -1;  // alternate signs
		double l1_alt = l0 - mstep_pen1_lagrange_cstr_func(l0, &md) /
			mstep_pen1_lagrange_cstr_derv(l0, &md);

		debug_msg(DEBUG_III, fxn_debug, "l0: %15.3f, l1: %15.3f; l1_alt: %15.3f;\n",l0, l1, l1_alt);

		if (l1 > 0)
			md.signs[arg_max_idx] = 1;
		//else if(l1_alt < lamb_lb){
		//   md.signs[arg_max_idx] = 1;
		//}
	} else {
		l0 = max_xi - md.rho + 1e-6;
	}

	// \lambda_0 = md.rho/md.omega
	double lambda = mstep_newton(_plfunc, _plderv, l0, 1e-6, 100, &md);

	if (max_xi < md.rho)
		debug_msg(DEBUG_I, fxn_debug, "lambda: %15.3f\n",lambda);


	//if(isnan(lambda)){
	//	debug_msg(DEBUG_I, fxn_debug, "Try again with another initialize value !\n");
	//	lambda = mstep_newton(_plfunc, _plderv, md.rho/md.omega, 1e-6, 100, &md);
	//}

	if (isnan(lambda)) { /* Use MLE instead */
		debug_msg(DEBUG_I, fxn_debug, "lambda for the %d th UMI: %15.3f\n", s, lambda);
		// return NAN;
		fprintf(stderr,"Fail to update lambda, use MLE instead of MPLE.\n");
		//debug_msg(DEBUG_I, fxn_debug, "ll penalty: %15.3f\n",lls);

		unsigned int idx;
		lls = 0.;
		for (int base = 0; base < md.n_xi; ++base) {
			idx = start+md.idx[base];
			x_s[idx]/= (sum_xi+ 1e-100);  // avoid /0
			lls+= md.rho * log1p(x_s[idx]/md.omega);
		}

		//for(unsigned int base = 0; base < md.n_xi; ++base)
		//		fprintf(stderr, "Est gamma: %15.3f", md.xi[base]);

		// fprintf(stderr, "\n");
		return lls;
	}

	double sum_tp = 0.0;
	lls = 0.0;
	for (int base = 0; base < md.n_xi; ++base) {
		/* compute the MPLE for transition prob. */
		double _xi = md.xi[base];
		double b = md.omega * lambda - _xi + rho;
		double sign = md.signs[base];
		double tp;

		tp = (-b + sign * sqrt(b * b + 4 * lambda * md.omega * _xi)) /
			 (2 * lambda);

		sum_tp += tp;
		lls += md.rho * log1p(tp/md.omega);
		x_s[start+md.idx[base]] = tp;
		if (isnan(tp) || tp < 0) {
			debug_msg(DEBUG_I, fxn_debug, "%d tp:%15.3f\n",md.idx[base],tp);
			debug_msg(DEBUG_I, fxn_debug, "lls :%15.3f\n",lls);
		}
		//if(s == 326)
		//	debug_msg(DEBUG_I, fxn_debug, "gamma :%15.3f; e: %15.3f\n", tp, _xi);
		if (max_xi < md.rho)
			debug_msg(DEBUG_III, fxn_debug, "_xi, %15.3f, est_x: %15.3f", _xi, tp);
	}
	if (isnan(lls))
		debug_msg(DEBUG_I, fxn_debug, "sum_tp:%15.3f\n",sum_tp);

	if (md.xi)
		free(md.xi);
	if (md.idx)
		free(md.idx);
	if (md.signs)
		free(md.signs);
	md.xi = NULL;
	md.idx = NULL;
	md.signs = NULL;

	return lls;
}/* MPLE_gamma */

int log_vector(double *x, unsigned int len)
{

	for (unsigned int i = 0; i < len; ++i) {
		if(x[i] > 0.)
			x[i] = log(x[i]);
		else
			x[i] = -INFINITY;
	}
	return NO_ERROR;
}/* log_vector */

int normalize(size_t n, unsigned int K, double * Emat)
{

	double max, sum;
	unsigned int idx;
	int fxn_debug = ABSOLUTE_SILENCE;

	for (unsigned int i = 0; i < n; ++i) {
		max = -INFINITY;
		for (unsigned int k = 0; k < K; k++) {
			idx = k * n + i;
			if (max < Emat[idx])
				max = Emat[idx];
		}
		if (i<10)
			debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, max);
		sum = 0;
		for (unsigned int k = 0; k < K; ++k) {
			idx = k * n + i;
			Emat[idx] = exp(Emat[idx] - max);
			sum += Emat[idx];
		}
		if (i<10)
			debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, sum);
		for (unsigned int k = 0; k < K; ++k) {
			idx = k * n + i;
			Emat[idx] /= sum;
			Emat[idx] = log(Emat[idx]);
		}
	}

	return NO_ERROR;
}/* normalize */


/* Assgin to the maximum maginal likelihood */
int reads_assign_sparse(model *mod, run_info *ri, size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI)
{


	int err = NO_ERROR;
	int fxn_debug = ABSOLUTE_SILENCE;
	double *rll = NULL;
	double *ull = NULL;

	rll = calloc(K, sizeof *rll);
	ull = calloc(K_UMI, sizeof *ull);
	if (!rll || !ull)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "reads assign: rll, ull");

	unsigned int idx;
	double rmax, umax;
	unsigned int rmax_id = 0, umax_id = 0;

	for (unsigned int k = 0; k < K; ++k)
		ri->optimal_cluster_size[k] = 0;

	for (unsigned int b = 0; b < K_UMI; ++b)
		ri->UMI_cluster_size[b] = 0;

	for (unsigned int i = 0; i < sample_size; ++i) {
		for (unsigned int j = 0; j < topN; ++j) {
			idx = i * topN + j;
			rll[mod->E2_sparse_hap_id[idx]] += mod->E2_sparse_value[idx];
			ull[mod->E2_sparse_umi_id[idx]] += mod->E2_sparse_value[idx];
		}

		rmax = 0.;
		for (unsigned int j = 0; j < K; ++j) {
			if (rll[j] > rmax) {
				rmax = rll[j];
				rmax_id = j;
			}
			rll[j] = 0;
		}

		ri->optimal_cluster_id[i] = rmax_id;
		ri->optimal_cluster_size[rmax_id]++;
		ri->optimal_cluster_ll[i] = rmax;

		if (i<10)
			debug_msg(DEBUG_I, fxn_debug, "%5d th read: %5d\n",i, rmax_id);


		umax = 0.;
		for (unsigned int j = 0; j < K_UMI; ++j) {
			if (ull[j] > umax) {
				umax = ull[j];
				umax_id = j;
			}
			ull[j] = 0;
		}

		ri->UMI_cluster_id[i] = umax_id;
		ri->UMI_cluster_size[umax_id]++;
		ri->UMI_cluster_ll[i] = umax;

		if (i<10)
			debug_msg(DEBUG_I, fxn_debug, "%5d th barcodes: %5d\n",i, umax_id);

	}

	/* free space */
	if (rll)
		free(rll);
	if (ull)
		free(ull);

	return err;
}/* reads_assign_sparse */

/* Assgin to the maximum joint likelihood */
int reads_assign_optimal(model *mod, run_info *ri, size_t sample_size, unsigned int topN,
					unsigned int K, unsigned int K_UMI)
{


	int err = NO_ERROR;

	unsigned int idx;

	for (unsigned int k = 0; k < K; ++k)
		ri->optimal_cluster_size[k] = 0;

	for (unsigned int b = 0; b < K_UMI; ++b)
		ri->UMI_cluster_size[b] = 0;

	for (unsigned int i = 0; i < sample_size; ++i) {
		idx = i * topN;
		ri->optimal_cluster_id[i] = mod->E2_sparse_hap_id[idx];
		ri->optimal_cluster_size[mod->E2_sparse_hap_id[idx]]++;

		ri->UMI_cluster_id[i] = mod->E2_sparse_umi_id[idx];
		ri->UMI_cluster_size[mod->E2_sparse_umi_id[idx]]++;

		ri->optimal_cluster_ll[i] = mod->E2_sparse_value[idx];

	}

	return err;
}/* reads_assign_optimal */


/* gam: gamma ; b : phi; xij = sum e_ij; thresh: rho */
double mstep_pen1_lagrange_cstr_func(double lambda, void *fdata)
{
	mstep_data *md = fdata;
	//mstep_data_t *md = fdata;
	double gam = md->omega;
	//int32_t sgns = md->signs;
	int n_xi = md->n_xi;
	//int sgn_off = n_xi > 4;

	register double sum_root = 0.0;
	for (int j = 0; j < n_xi; j++) {
		double thresh = md->rho;
		double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) {
			fprintf(stderr, "Delta < 0 in func\n");
			return NAN;
		}

		delta = sqrt(delta);

		/* try mod->signs later */
		double sign = md->signs[j]; //(double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		sum_root += (-b + sign * delta)/(2 * lambda);
	}

	return sum_root - 1;
}


double mstep_pen1_lagrange_cstr_derv(double lambda, void *fdata)
{
	//register mstep_data_t *md = fdata;
	mstep_data *md = fdata;
	register double gam = md->omega;
	//register int32_t sgns = md->signs;
	int n_xi = n_xi;
	//int sgn_off = n_xi > 4;

	register double derv = 0.0;
	for (int j = 0; j < md->n_xi; j++) {
		register double thresh =  md->rho; //eta / log1p(1.0/gam);
		register double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0) {
			fprintf(stderr, "Delta < 0 in derv\n");
			return NAN;
		}

		delta = sqrt(delta);
		double sign = md->signs[j]; //(double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
		derv += (lambda * (-gam + sign / delta * gam * (b + 2 * xij))) -
			(-b + sign * delta);
	}

	return derv / (2 * lambda * lambda);
}

/* ----- Root finding algorithms ----- */
double mstep_newton(double (*fx)(double x, void *data),
		double (*fderv)(double x, void *data), double x0, double eps,
		int maxiter, void *fdata)
{
	int iter = 1;
	double x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

	if (isnan(x1))
		fprintf(stderr, "%5d iter: x0: %15.3f,x1:%15.3f\n", iter, x0, x1);

	if (isnan(x1)) return NAN;

	while (fabs(x1 - x0) > eps && iter < maxiter) {

		x0 = x1;
		x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

		if (isnan(x1))
			fprintf(stderr, "%5d iter: x0: %15.3f,x1:%15.3f\n", iter, x0, x1);

		if (isnan(x1))
			return NAN;

		iter++;
	}

	return x1;
}

double mstep_pen1_lambda_support(void *fdata)
{
	mstep_data *md = fdata;
	register double _omega = md->omega;
	register int n_xi = md->n_xi;
	register double _rho = md->rho;

	register double max_lamb_lb = -INFINITY;
	for (int j = 0; j < n_xi; j++) {
		register double xi_j = md->xi[j];
		register double _lb = (-(xi_j + _rho) + 2 * sqrt(xi_j * _rho)) / _omega;
		if (_lb > max_lamb_lb)
			max_lamb_lb = _lb;
	}

	return max_lamb_lb;
}

int trans_expect_UMIs(options *opt, data *dat, data_t *seeds_UMI,
	double *error_profile, double *trans_prob, int ends_free)
{

	int err = NO_ERROR;
	hash *s;
	unsigned int u = 0;
	double eik;
	double l1third = 1.0/3;


	if (!opt->band || opt->nw_align == NO_ALIGNMENT) {

		for (unsigned int r = 0; r < dat->sample_size; ++r) {

			for (unsigned int b = 0; b < opt->K_UMI; b++) {

				eik = 0.;
				unsigned char * hap_seq = &seeds_UMI[b *opt->UMI_length];

				for (unsigned int j = 0; j < opt->UMI_length; j++) {

					if (error_profile) {
						if (opt->err_encoding == STD_ENCODING)  // make sure opt->err_encoding is the same as mod->err_encoding
							eik += translate_error_STD_to_XY(
								error_profile,
								dat->n_quality, hap_seq[j],
								dat->dmatU[r][j],
								dat->qmatU[r][j]);
						else if (opt->err_encoding == XY_ENCODING)
							eik += error_profile[(NUM_NUCLEOTIDES
								* hap_seq[j] + dat->dmatU[r][j])
								* dat->n_quality
								+ dat->qmatU[r][j]];
					} else {
						//double ep = adj * error_prob(dat->fdata, dat->qmat[r][j]);
						double ep = dat->error_prob[dat->qmatU[r][j]];
						if (dat->dmatU[r][j] == hap_seq[j] )
							eik += log(1 - ep);
						else
							eik += log(ep) + l1third;
					}
				}


				trans_prob[b * dat->sample_size + r] = eik;

			}

		}
		return err;
	}

	double adj_trunpois = ppois(dat->max_read_length, dat->max_read_length * opt->indel_error, 1, 1);

	for (s = dat->UMI_count; s != NULL; s = s->hh.next) {

		if (u == dat->hash_UMI_length)
			return mmessage(ERROR_MSG, INTERNAL_ERROR, "exceed the length");

		unsigned char *read = s->sequence;
		unsigned int rlen = opt->UMI_length;

		unsigned int count = s->count; // num. of reads wih unique seq
		size_t *idx_array = s->idx_array;  // idx of reads

		for (unsigned int b = 0; b < opt->K_UMI; ++b) {

			unsigned char *hap_seq = &seeds_UMI[b *opt->UMI_length];
			size_t alen = opt->UMI_length;
			unsigned int nindels = 0;
			unsigned int nmismatch = 0;

			unsigned char **aln = nwalign(hap_seq, read,
				(size_t) opt->UMI_length, (size_t) rlen,
				opt->score, opt->gap_p, opt->band, 1, NULL,
								&err, &alen);

			/* count for number of indels */
			ana_alignment(aln, alen, rlen, &nindels,
					&nmismatch, opt->ends_free, opt->info); // need further check

			for (unsigned int r = 0; r<count;++r) {

				/*
				if (idx_array[r] < 10 && b == 0 ) {
					for (size_t j = 0; j < alen; ++j) {
						fprintf(stderr, "%c", aln[0][j] == '-'
							? '-' : xy_to_char[(int) aln[0][j]]);
					}
					fprintf(stderr, "\n");
					for (size_t j = 0; j < alen; ++j) {
						fprintf(stderr, "%c", aln[1][j] == '-'
							? '-' : xy_to_char[(int) aln[1][j]]);
					}
					fprintf(stderr, "\n");
				}
				*/

				trans_prob[b * dat->sample_size + idx_array[r]] = trans_nw(opt, aln,
					alen, nmismatch, nindels, error_profile, opt->err_encoding,
					dat->qmatU[idx_array[r]], dat->n_quality, adj_trunpois,
								rlen, dat->error_prob, ends_free);

			}

			if (aln) {
				free(aln[0]);
				free(aln[1]);
				free(aln);
				aln = NULL;
			}
			if (err)
				return err;
		}

		u++;
	}

	return err;
}
