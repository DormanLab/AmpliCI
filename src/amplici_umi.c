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
double MPLE_gamma_s(double *x_s, unsigned int K, int *err, unsigned int s,double pho, double omega);
int trans_expect_UMIs(options *opt, data *dat, data_t *seeds_UMI, double *error_profile, double *trans_prob);


/* main model for UMI */
int EM_algorithm(options *opt, data *dat, model *mod, initializer *ini, run_info *ri){
   int err = NO_ERROR;
   double rdelta;
   int EM_iter;

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

    do{

    mmessage(INFO_MSG, NO_ERROR, "The E step of the %5d iteration...\n", EM_iter);

   log_vector(mod->gamma, opt->K*opt->K_UMI);
   log_vector(mod->eta, opt->K_UMI);
    
    mod->ll_UMI = E_step(mod, dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI, &err);
    mod->ll_UMI -= mod->penalty_ll;
    if(err)
        return err;


    rdelta = (mod->pll_UMI - mod->ll_UMI)/mod->pll_UMI;
    mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (%15.3f - %15.3f, %5.3e)\n", EM_iter, 
                    mod->ll_UMI,mod->ll_UMI+mod->penalty_ll, mod->penalty_ll, rdelta);
   

    mmessage(INFO_MSG, NO_ERROR, "The M step of the %5d iteration...\n", EM_iter);

    /* update gamma and eta */
    if((err = M_step(opt,mod,dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI)))
        return err;


    mod->pll_UMI = mod->ll_UMI;
    EM_iter ++;

    if (rdelta < opt->epsilon && rdelta > 0)
		break;


    if (EM_iter > opt->n_iter_amplici){
		mmessage(WARNING_MSG, INTERNAL_ERROR,
			"exceed max interations,%i\n", EM_iter);
        break;
    }

    }while(1);

    mmessage(INFO_MSG, NO_ERROR, "Final reads assignment.....\n");

    /* reads assignment */
    /*
    if((err = reads_assign_sparse(mod, ri, dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI)))
        return err;
    */

    /* assignment with the combination with the joint maximum log likelihood */
    reads_assign_optimal(mod, ri, dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI);
    

     /* print and check */
    FILE *fp = NULL;
    fp = fopen(opt->outfile_base, "w");
    if (!fp){
        return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
                        opt->outfile_base);        
    }

    fprintf(fp, "log likelihood: %f\n", mod->ll_UMI);

    fprintf(fp, "K: %i\n", opt->K);

    fprintf(fp, "assignments: ");
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

    fprintf(fp,"reads posterior likelihood: ");
	fprint_doubles(fp, ri->optimal_cluster_ll, dat->sample_size, 3, 1);

    fprintf(fp,"Eta: ");
	fprint_doubles(fp, mod->eta, opt->K_UMI, 6, 1);

    fprintf(fp,"Gamma: ");
    fprint_vectorized_matrix(fp, mod->gamma,
				opt->K_UMI, opt->K, 1);
    
    fclose(fp);

    mmessage(INFO_MSG, NO_ERROR, "Output the final result file: "
						"%s \n", opt->outfile_base);

    return err;
}/* EM_algorithm */


/* set initial value for gamma and eta */
double ini_eta_gamma(options *opt, initializer *ini, double *gamma, double *eta, size_t sample_size)
{

    //int err = NO_ERROR;
    unsigned int K = opt->K;
    int fxn_debug = ABSOLUTE_SILENCE;
    double llp_sum;
    double rho = 1.01;
    double omega = 1e-15;
    double epsilon = 1e-40;


    for (size_t i = 0; i < sample_size; i++)
    {
        if (ini->reads_umi_id[i] > -1)
        {
            eta[ini->reads_umi_id[i]] += 1.0;
            if (ini->reads_hap_id[i] > -1)
            {
                gamma[ini->reads_umi_id[i] * K + ini->reads_hap_id[i]] += 1.0;
            }
        }
    }

    double eta_sum = 0.;
    unsigned int idx;

    llp_sum = 0.;
    for (unsigned int b = 0; b < opt->K_UMI; ++b)
    {
        eta_sum += eta[b];

        debug_msg(DEBUG_II, fxn_debug, "%15.3f,",eta[b]);
        debug_msg(DEBUG_II, fxn_debug, "\n");

        double gamma_sumh = 0.;
        for (unsigned int k = 0; k < K; ++k)
            gamma_sumh += gamma[b * K + k];

        for (unsigned int k = 0; k < K; ++k){
            idx = b * K + k;
            debug_msg(DEBUG_II, fxn_debug, "%15.3f / %15.3f\n,",gamma[idx],gamma_sumh);
            gamma[idx] /=  (gamma_sumh+ epsilon); // avoid /0
            debug_msg(DEBUG_III, fxn_debug, "%15.3f,",gamma[idx]);  
            double llp = rho * log1p(gamma[idx]/omega);
            llp_sum += llp;
            if(isnan(llp))
                fprintf(stderr, "(%d, %d) : gamma = %15.7f, ll = %15.7f \n", b, k, gamma[idx], llp);
        }
          debug_msg(DEBUG_III, fxn_debug, "\n");  
    }

    for (unsigned int b = 0; b < opt->K_UMI; ++b){
        eta[b] /= eta_sum;
    }

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
    if (opt->use_error_profile && mod->error_profile)
    {
        error_profile = mod->error_profile;
        debug_msg(DEBUG_II, fxn_debug, "Use error profile. \n");
    }

     mmessage(INFO_MSG, NO_ERROR, "Calculating the trans prob of reads......\n");


    /* For acceleration purpose. Do not allow large gap in amplification and sequencing */
    opt->band = 5;
   
    /* transition probability for reads */
    if ((err = trans_expectation(opt, dat, ini, error_profile,
                                 mod->adj_trunpois, mod->eik)))
        return err;


    /* transition probability for UMIs */
    /* implement a naive method, waiting for an efficient one */

    mmessage(INFO_MSG, NO_ERROR, "Calculating the trans prob of UMIs......\n");

     /* [TODO] Need further investigation. It is good to set high gap 
            penalty to reduce gap when aligning UMIs */
    opt->band = 2;   // for UMIs 
    opt->gap_p = -20;
    
    /* transition probability for UMIs */
    if((err = trans_expect_UMIs(opt, dat, ini->seeds_UMI, error_profile, mod->eik_umi)))
        return err;
        
   
   /*
    normalize(dat->sample_size, opt->K, mod->eik);
    normalize(dat->sample_size, opt->K_UMI, mod->eik_umi);
    */
    

    return err;
}/* trans_hap_and_umi */

/* update E2_sparse */
double E_step(model *mod, size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI, int *err){

    double ll = 0.;
    unsigned int idx;
    unsigned int s;
    double max, cur_value;
    unsigned int max_idx;
    double maximum;
    double sum;
    int fxn_debug = ABSOLUTE_SILENCE;
    double read_ll;
    // double log_epsilon = - 40;
    
    /* [TODO] temp is sparse, think about a more efficient way to do */
    double *temp = NULL;
    temp = calloc(K* K_UMI, sizeof *temp);
    if(!temp){
        *err = MEMORY_ALLOCATION;
        mmessage(ERROR_MSG, MEMORY_ALLOCATION, "temp");
        return -INFINITY;
    }

    double *gamma = mod->gamma;
    double *eta = mod->eta;


    /* Find top N log joint probabilities for each read */
    for (unsigned int i = 0; i < sample_size; ++i){
        idx = i * topN;
        s = 0;
        
        /* Calculation and find the maximum value */
        max = -INFINITY;
        for (unsigned int k = 0; k < K; ++k){
            for (unsigned int b = 0; b < K_UMI; ++b){
                    if(gamma[b*K+k] > -INFINITY && eta[b] > -INFINITY){
                        cur_value = gamma[b*K + k] + eta[b] 
                            + mod->eik_umi[b*sample_size +i] + mod->eik[k*sample_size + i];
                   }else{
                        debug_msg(DEBUG_I, fxn_debug, "gamma == 0\n");
                        cur_value = -INFINITY;
                    }
                    temp[b*K+k] = cur_value;
                    if(cur_value > max){
                        max = cur_value;
                        mod->E2_sparse_value[idx + s] = cur_value;
                        mod->E2_sparse_hap_id[idx + s] = k;
                        mod->E2_sparse_umi_id[idx + s] = b;               
                    }
            }
        }
        if(i < 10){
         debug_msg(DEBUG_I, fxn_debug, "maximum value for the %5d reads: %15.3f\n",i, max);
         debug_msg(DEBUG_I, fxn_debug, "gamma : %15.3f\n",gamma[mod->E2_sparse_umi_id[idx + s]*K+mod->E2_sparse_hap_id[idx + s]]);
         debug_msg(DEBUG_I, fxn_debug, "eta : %15.3f\n",eta[mod->E2_sparse_umi_id[idx + s]]);
        debug_msg(DEBUG_I, fxn_debug, "EUMI : %15.3f\n",mod->eik_umi[mod->E2_sparse_umi_id[idx + s]*sample_size +i]);
        debug_msg(DEBUG_I, fxn_debug, "E : %15.3f\n",mod->eik[mod->E2_sparse_hap_id[idx + s]*sample_size +i]);
        }
        
        temp[mod->E2_sparse_umi_id[idx + s]*K+mod->E2_sparse_hap_id[idx + s]] = -INFINITY;
        s += 1;
        maximum = max;  // maximum value among top N
       
        while (s < topN)
        {    
            max = -INFINITY;
            for (unsigned int k = 0; k < K; ++k)
            {
                for (unsigned int b = 0; b < K_UMI; ++b)
                {   
                    cur_value = temp[b * K + k];
                    if (cur_value > max){
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
        for (unsigned int t = 0; t < topN; ++t){
            mod->E2_sparse_value[idx + t] = exp(mod->E2_sparse_value[idx + t] - maximum);
            sum += mod->E2_sparse_value[idx + t];
        }

        read_ll = log(sum) +maximum;
        ll += read_ll; 
        if(i < 10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, read_ll);

        for (unsigned int t = 0; t < topN; ++t)
            mod->E2_sparse_value[idx + t] /= sum;

    }
    debug_msg(DEBUG_III, fxn_debug, "\n");

    if(temp) free(temp);

    return ll;
}/* E_step */



/* estimate new eta and new gamma */
int M_step(options *opt, model *mod,size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI){
    int err = NO_ERROR;
    int fxn_debug = DEBUG_I;
    unsigned int idx; 
    double epsilon = 1e-40;
    int algo = opt->trans_penalty;

    double *eta =  mod->eta;
    double *gamma = mod->gamma;

   
    for (unsigned int b = 0; b < K_UMI; ++b){  
            eta[b] = 0.;
         for (unsigned int k = 0; k < K; ++k)
             gamma[b * K + k] = 0.;
         
    } 

    for (unsigned int i= 0; i< sample_size; ++i){
        for (unsigned int t = 0; t < topN;++t){
            idx = i *topN + t;
            gamma[mod->E2_sparse_umi_id[idx]* K + mod->E2_sparse_hap_id[idx]] += mod->E2_sparse_value[idx];
        }
    }

    for (unsigned int b = 0; b < K_UMI; ++b){  
         for (unsigned int k = 0; k < K; ++k)
             eta[b] += gamma[b * K + k]; 
    } 

    //for(unsigned int k = 0; k < K; ++k)
     //   fprintf(stderr, "Est gamma: %15.3f", gamma[k]);
    
   // debug_msg(DEBUG_I, fxn_debug, "\n");

    /* gamma */
    mod->penalty_ll = 0.;
    for (unsigned int b = 0; b < K_UMI; b++){

        if(algo == MPLE){
            /* for each row of gamma */
            double lls = MPLE_gamma_s(gamma, K, &err, b,opt->pho, opt->omega);
            if(err) return err;

            
            if (isnan(lls)){
                debug_msg(DEBUG_I, fxn_debug, "ll penalty: %15.3f\n",lls);
                
                // compute MLE instead [XY: is it fine ? ]
                for (unsigned int k = 0; k < K; k++){
                    idx = b * K + k;
                    gamma[idx] /= (eta[b] + epsilon); // avoid /0
                    double llp = opt->pho * log1p(gamma[idx]/opt->omega);
                    mod->penalty_ll += llp;
                }
                
            }else{
                mod->penalty_ll += lls;
            }
            
        }else{
            for (unsigned int k = 0; k < K; k++){
                idx = b * K + k;
                gamma[idx] /= (eta[b] + epsilon); // avoid /0
            }
        }
    }

   
    //for(unsigned int k = 0; k < K; ++k){
    //    fprintf(stderr, "Est gamma: %15.3f", gamma[k]);
    //}
    //debug_msg(DEBUG_I, fxn_debug, "\n");


    /* eta */
    for (unsigned int b = 0; b < K_UMI; b++){
        eta[b] /= sample_size;
    }

    return err;
}/* M_step */

/* Caculate the MPLE for gamma for each s and return the panalty turn */
double MPLE_gamma_s(double *x_s, unsigned int K, int *err, unsigned int s, double pho, double omega){

    int fxn_debug =DEBUG_I;
    double lls;
    int start = s*K;

    /* set up */
    mstep_data md = {
		.omega =omega,
        .pho = pho,
		.xi = NULL,
        .idx = NULL,
		.signs = 0,
		.n_xi = 0,
	};


    /* count the no-zeros number */
    for(unsigned int k =0 ; k < K; ++k)
        if(x_s[start+k] > 0)
            md.n_xi+= 1;

    if(md.n_xi == 0)
        return 0;

    if(md.n_xi == 1)
        return md.pho * log1p(1.0/md.omega);

    md.xi = calloc(md.n_xi, sizeof (*x_s));
    md.idx = calloc(md.n_xi, sizeof (*md.idx));

    if(!md.xi || !md.idx){
        *err = MEMORY_ALLOCATION;
        mmessage(ERROR_MSG, *err,"md\n");
        return 0;
    }
        
    unsigned int cnt = 0;
    for(unsigned int k = 0; k < K; ++k){
       if(x_s[start+k] > 0){
           md.xi[cnt] = x_s[start+k];
           md.idx[cnt] = k;
            cnt ++;
       }
    }

    msp_lcstr_func _plfunc = mstep_pen1_lagrange_cstr_func;
	msp_lcstr_derv _plderv = mstep_pen1_lagrange_cstr_derv;

    // \lambda_0 = md.pho/md.omega
    double lambda = mstep_newton(_plfunc, _plderv, 1, 1e-6, 100, &md);
    
    if (isnan(lambda)){
        debug_msg(DEBUG_I, fxn_debug, "lambda for the %d th UMI: %15.3f\n", s, lambda);
        return NAN;
    }

    double sum_tp = 0.0;
    lls = 0.0;
    for (unsigned int base = 0; base < md.n_xi; ++base)
    {
        double rho = md.pho;
        /* compute the MPLE for transition prob. */
        double _xi = md.xi[base];
        double b = md.omega * lambda - _xi + rho;
        double sign = 1;
        double tp;

        tp = (-b + sign *
                       sqrt(b * b + 4 * lambda * md.omega * _xi)) /
             (2 * lambda);

        sum_tp += tp;
        lls += rho * log1p(tp/md.omega);
        x_s[start+md.idx[base]] = tp;
        if (isnan(tp)){
            debug_msg(DEBUG_I, fxn_debug, "%d tp:%15.3f\n",md.idx[base],tp);
            debug_msg(DEBUG_I, fxn_debug, "lls :%15.3f\n",lls);
        }
    }
    if (isnan(lls))
        debug_msg(DEBUG_I, fxn_debug, "sum_tp:%15.3f\n",sum_tp);  

    if(md.xi) free(md.xi);
    if(md.idx) free(md.idx);
    md.xi = NULL;
    md.idx = NULL;

    return lls;
}/* MPLE_gamma */

int log_vector(double *x, unsigned int len){

    for(unsigned int i = 0; i < len; ++i){
        if(x[i] > 0.){
            x[i] = log(x[i]);
        }else{
            x[i] = -INFINITY;
        }
    }
    return NO_ERROR;
}/* log_vector */

int normalize(size_t n, unsigned int K, double * Emat){

    double max, sum;
    unsigned int idx;
    int fxn_debug = ABSOLUTE_SILENCE;

     for (unsigned int i = 0; i < n; ++i){
        max = -INFINITY;
        for (unsigned int k = 0; k < K; k++) {
			idx = k * n + i;
			if (max < Emat[idx])
				max = Emat[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, max);
		sum = 0;
		for (unsigned int k = 0; k < K; ++k) {
			idx = k * n + i;
			Emat[idx] = exp(Emat[idx] - max);
			sum += Emat[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, sum);
        for (unsigned int k = 0; k < K; ++k){
            idx = k * n + i;
			Emat[idx] /= sum;
            Emat[idx] = log(Emat[idx]);
        }
    }

    return NO_ERROR;
}/* normalize */


/* Assgin to the maximum maginal likelihood */
int reads_assign_sparse(model *mod, run_info *ri, size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI){


    int err = NO_ERROR;
    int fxn_debug = ABSOLUTE_SILENCE;
    double *rll = NULL;
    double *ull = NULL;

    rll = calloc(K, sizeof *rll);
    ull = calloc(K_UMI, sizeof *ull);
    if(!rll || !ull)
        return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "reads assign: rll, ull");
                   
    unsigned int idx;
    double rmax, umax;
    unsigned int rmax_id, umax_id;

    for (unsigned int k = 0; k < K; ++k)
        ri->optimal_cluster_size[k] = 0;

    for (unsigned int b = 0; b < K_UMI; ++b)
        ri->UMI_cluster_size[b] = 0;

    for (unsigned int i = 0; i < sample_size; ++i)
    {
        for (unsigned int j = 0; j < topN; ++j)
        {
            idx = i * topN + j;
            rll[mod->E2_sparse_hap_id[idx]] += mod->E2_sparse_value[idx];
            ull[mod->E2_sparse_umi_id[idx]] += mod->E2_sparse_value[idx];
        }

        rmax = 0.;
        for (unsigned int j = 0; j < K; ++j)
        {
            if (rll[j] > rmax)
            {
                rmax = rll[j];
                rmax_id = j;
            }
            rll[j] = 0;
        }

        ri->optimal_cluster_id[i] = rmax_id;
        ri->optimal_cluster_size[rmax_id]++;
        ri->optimal_cluster_ll[i] = rmax;

        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %5d\n",i, rmax_id);


        umax = 0.;
        for (unsigned int j = 0; j < K_UMI; ++j)
        {
            if (ull[j] > umax)
            {
                umax = ull[j];
                umax_id = j;
            }
            ull[j] = 0;
        }

        ri->UMI_cluster_id[i] = umax_id;
        ri->UMI_cluster_size[umax_id]++;
        ri->UMI_cluster_ll[i] = umax;

        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th barcodes: %5d\n",i, umax_id);

    }

    /* free space */
    if(rll)free(rll);
    if(ull)free(ull);

    return err;
}/* reads_assign_sparse */

/* Assgin to the maximum joint likelihood */
int reads_assign_optimal(model *mod, run_info *ri, size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI){


    int err = NO_ERROR;

    unsigned int idx;
   
    for (unsigned int k = 0; k < K; ++k)
        ri->optimal_cluster_size[k] = 0;

    for (unsigned int b = 0; b < K_UMI; ++b)
        ri->UMI_cluster_size[b] = 0;

    for (unsigned int i = 0; i < sample_size; ++i)
    {
        idx = i * topN;
        ri->optimal_cluster_id[i] = mod->E2_sparse_hap_id[idx];
        ri->optimal_cluster_size[mod->E2_sparse_hap_id[idx]]++;

        ri->UMI_cluster_id[i] = mod->E2_sparse_umi_id[idx];
        ri->UMI_cluster_size[mod->E2_sparse_umi_id[idx]]++;

        ri->optimal_cluster_ll[i] = mod->E2_sparse_value[idx];

    }

    return err;
}/* reads_assign_sparse */


/* gam: gamma ; b : phi; xij = sum e_ij; thresh: pho */
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
		double thresh = md->pho;
		double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0){ 
            fprintf(stderr, "Delta < 0 in func\n");
            return NAN;
        }

		delta = sqrt(delta);

        /* try mod->signs later */
		double sign = 1; //(double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
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
		register double thresh =  md->pho; //eta / log1p(1.0/gam);
		register double xij = md->xi[j];

		double b = gam * lambda - xij + thresh;
		double mc = gam * xij;

		double delta = b * b + 4 * lambda * mc;
		/* sanity check */
		if (delta < 0){ 
            fprintf(stderr, "Delta < 0 in derv\n");
            return NAN;
        }

		delta = sqrt(delta);
		double sign = 1; //(double) (((((sgns >> j) & 1) + sgn_off) << 1) - 1);
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
        if (isnan(x1))
            fprintf(stderr, "%5d iter: x0: %15.3f,x1:%15.3f\n", iter, x0, x1);

		x0 = x1;
		x1 = x0 - fx(x0, fdata) / fderv(x0, fdata);

		if (isnan(x1)) return NAN;

		iter++;
	}

	return x1;
}

int trans_expect_UMIs(options *opt, data *dat, data_t *seeds_UMI, double *error_profile, double *trans_prob){

    int err = NO_ERROR;
    hash *s;
	unsigned int u = 0;

    double adj_trunpois = ppois(dat->max_read_length, dat->max_read_length * opt->indel_error, 1, 1);

	for (s = dat->UMI_count; s != NULL; s = s->hh.next) {

		if (u == dat->hash_UMI_length)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
							"exceed the length");

        unsigned char *read = s->sequence;
		unsigned int rlen = opt->UMI_length;

		unsigned int count = s->count; // num. of reads wih unique seq
		size_t *idx_array = s->idx_array;  // idx of reads

       for (unsigned int b = 0; b < opt->K_UMI; ++b){
           
           unsigned char *hap_seq = &seeds_UMI[b *opt->UMI_length];

            size_t alen = opt->UMI_length;
		    unsigned int nindels = 0;
			unsigned int nmismatch = 0;

				unsigned char **aln = nwalign(hap_seq, read,
				(size_t) opt->UMI_length,
				(size_t) rlen,
				opt->score, opt->gap_p, opt->band, 1, NULL,
									&err, &alen);

				/* count for number of indels */
				ana_alignment(aln, alen, rlen, &nindels, 
						&nmismatch, opt->info); // need further check 

				for(unsigned int r = 0; r<count;++r){

					trans_prob[b * dat->sample_size + idx_array[r]] = trans_nw(opt, aln,
						alen, nmismatch, nindels, error_profile, opt->err_encoding,
						dat->qmatU[idx_array[r]], dat->n_quality, adj_trunpois,
										rlen, dat->error_prob);
				
				}

                if(aln){
					free(aln[0]);
					free(aln[1]);
					free(aln);
					aln = NULL;
				}
				if(err)
					return err;
       }

		u++;
	}

    return err;
}
