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

#define MATHLIB_STANDALONE 1
#include <Rmath.h>


int trans_hap_and_umi(options *opt, data *dat, initializer *ini, model *mod);
double E_step(model *mod, size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI, int *err);
int M_step(model *mod,size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI);
int reads_assign_sparse(model *mod, run_info *ri, size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI);

/* main model for UMI */
int EM_algorithm(options *opt, data *dat, model *mod, initializer *ini, run_info *ri){
   int err = NO_ERROR;
   double rdelta;
   int EM_iter;

   /* calculate the transition prob */
   if ((err = trans_hap_and_umi(opt, dat, ini, mod)))  
        return err;


    /* initialize gamma and eta */
    ini_eta_gamma(opt, ini, mod->gamma, mod->eta, dat->sample_size);
     

    /* when ll relative diff > 1e-6 */ 

    /* simple E step, output E_id and E_value and ll given pgamma and peta */

    EM_iter = 0;

    do{

    mmessage(INFO_MSG, NO_ERROR, "The E step of the %5d iteration...\n", EM_iter);


    /* rewrite below as a function */
    /* 
    unsigned int idx;
     for (unsigned int b = 0; b < opt->K_UMI; b++){
        for (unsigned int k = 0; k < opt->K; k++){
             idx = b * opt->K + k;
             if(mod->gamma[idx]){
                mod->gamma[idx] = log(mod->gamma[idx]);   // store in log version
            }else{
                mod->gamma[idx] = -INFINITY;   
            }
        }
         if(mod->eta[b]){
            mod->eta[b] = log(mod->eta[b]);   // store in log version 
        }else{
            mod->eta[b] = -INFINITY;
        }
    }
    */
   log_vector(mod->gamma, opt->K*opt->K_UMI);
   log_vector(mod->eta, opt->K_UMI);
    
    mod->ll_UMI = E_step(mod, dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI, &err);
    if(err)
        return err;


    rdelta = (mod->pll_UMI - mod->ll_UMI)/mod->pll_UMI;
    mmessage(INFO_MSG, NO_ERROR, "%5d: %15.3f (%5.3e)\n", EM_iter, mod->ll_UMI, rdelta);
   

    mmessage(INFO_MSG, NO_ERROR, "The M step of the %5d iteration...\n", EM_iter);

    /* update gamma and eta */
    if((err = M_step(mod,dat->sample_size, opt->topN, 
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
    if((err = reads_assign_sparse(mod, ri, dat->sample_size, opt->topN, 
                    opt->K, opt->K_UMI)))
        return err;

     /* print and check */
    FILE *fp = NULL;
    fp = fopen(opt->outfile_base, "w");
    if (!fp){
        return mmessage(ERROR_MSG, FILE_OPEN_ERROR,
                        opt->outfile_base);        
    }

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
int ini_eta_gamma(options *opt, initializer *ini, double *gamma, double *eta, size_t sample_size)
{

    int err = NO_ERROR;
    unsigned int K = opt->K;
    int fxn_debug = ABSOLUTE_SILENCE;

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
            gamma[idx] /=  gamma_sumh;
            debug_msg(DEBUG_III, fxn_debug, "%15.3f,",gamma[idx]);  
        }
          debug_msg(DEBUG_III, fxn_debug, "\n");  
    }

    for (unsigned int b = 0; b < opt->K_UMI; ++b){
        eta[b] /= eta_sum;
    }
        

    return err;
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


   
    /* [TODO] Need further investigation. It is good to set high gap penalty to reduce gap in alignment */
    opt->gap_p = -20;

    /* For acceleration purpose. Do not allow large gap in amplification and sequencing */
    opt->band = 5;
   
    /* transition probability for reads */
    if ((err = trans_expectation(opt, dat, ini, error_profile,
                                 mod->adj_trunpois, mod->eik)))
        return err;


    /* transition probability for UMIs */
    /* implement a naive method, waiting for an efficient one */

     mmessage(INFO_MSG, NO_ERROR, "Calculating the trans prob of UMIs......\n");
    
    size_t alen;
    unsigned char **aln;
    unsigned int nindels,nmismatch;
    double adj_trunpois = ppois(dat->max_read_length, dat->max_read_length * opt->indel_error, 1, 1);

    
    /* seeks for a more efficient way */
    for (unsigned int r = 0; r < dat->sample_size; ++r)
    {

        for (unsigned int b = 0; b < opt->K_UMI; b++)
        {


            /* always find problems when aligning barcodes. maybe just use alignmnet free method ? */

            alen = opt->UMI_length;
            nindels = 0;
			nmismatch = 0;
            aln = nwalign(&ini->seeds_UMI[b * opt->UMI_length], dat->dmatU[r],
                                          (size_t)opt->UMI_length,
                                          (size_t)opt->UMI_length,
                                          opt->score, opt->gap_p, 2, 1, NULL,
                                          &err, &alen);        // set band = 2 


            ana_alignment(aln,alen, opt->UMI_length, &nindels, 
						&nmismatch, opt->info);

           /* 
             if(r ==0 && b < 10 ){
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
            
            

            mod->eik_umi[b * dat->sample_size + r] = trans_nw(opt, aln,
                                                              alen,nmismatch, nindels, error_profile, opt->err_encoding,
                                                              dat->qmatU[r], dat->n_quality, adj_trunpois,
                                                              opt->UMI_length, dat->error_prob);
            
            //if(r == 0 && b < 10 ){
                 //   for (size_t j = 0; j < alen; ++j){
                  //      fprintf(stderr, " %8.2e\n", error_profile[(
					//	NUM_NUCLEOTIDES
				//		* aln[0][j] + aln[1][j])
			//			* dat->n_quality + dat->qmatU[r][j]]);
              //          fprintf(stderr, " %c\n", dat->qmatU[r][j]+dat->fdata->min_quality);
                //    }
                //    fprintf(stderr, "\n");
             //       fprintf(stderr, " %8.2e\n", mod->eik_umi[b * dat->sample_size + r]);
			//}
            
        
            /* free */
				if(aln){
					free(aln[0]);
					free(aln[1]);
					free(aln);
					aln = NULL;
				}
				if(err)
					return err;
        }
    }
    
    /*
    if(opt->trans_matrix){
		FILE *fp = fopen(opt->trans_matrix, "w");
		if (!fp)
			return mmessage(ERROR_MSG, FILE_OPEN_ERROR, opt->trans_matrix);
		//fprint_vectorized_matrix(fp,mod->eik,dat->sample_size,opt->K,0);  not work 
		for (size_t i = 0; i < dat->sample_size; ++i) {
		    fprintf(fp, "%3lu", i);
			for (unsigned int j = 0; j < opt->K_UMI; ++j) {
				fprintf(fp, " %8.2e", mod->eik_umi[j*dat->sample_size + i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
    */

    normalize(dat->sample_size, opt->K, mod->eik);
    normalize(dat->sample_size, opt->K_UMI, mod->eik_umi);

     /* normalize and restored in log version */
     /*
    double max, sum;
    unsigned int idx;
    for (unsigned int i = 0; i < dat->sample_size; i++){
        max = -INFINITY;
        for (unsigned int k = 0; k < opt->K; k++) {
			idx = k * dat->sample_size + i;
			if (max < mod->eik[idx])
				max = mod->eik[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, max);
		sum = 0;
		for (unsigned int k = 0; k < opt->K; ++k) {
			idx = k * dat->sample_size + i;
			mod->eik[idx] = exp(mod->eik[idx] - max);
			sum += mod->eik[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th read: %15.3f\n",i, sum);
        for (unsigned int k = 0; k < opt->K; ++k){
            idx = k * dat->sample_size + i;
			mod->eik[idx] /= sum;
            mod->eik[idx] = log(mod->eik[idx]);
        }
    }
    */

    /* UMIs */
    /*
    for (unsigned int i = 0; i < dat->sample_size; i++){
        max = -INFINITY;
        for (unsigned int k = 0; k < opt->K_UMI; k++) {
			idx = k * dat->sample_size + i;
			if (max < mod->eik_umi[idx])
				max = mod->eik_umi[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th barcodes: %15.3f\n",i, max);
		sum = 0;
		for (unsigned int k = 0; k < opt->K_UMI; ++k) {
			idx = k * dat->sample_size + i;
			mod->eik_umi[idx] = exp(mod->eik_umi[idx] - max);
			sum += mod->eik_umi[idx];
		}
        if(i<10)
            debug_msg(DEBUG_I, fxn_debug, "%5d th barcodes: %15.3f\n",i, sum);
        for (unsigned int k = 0; k < opt->K_UMI; ++k){
            idx = k * dat->sample_size + i;
			mod->eik_umi[idx] /= sum;
            mod->eik_umi[idx] = log(mod->eik_umi[idx]);
        }
    }
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
int M_step(model *mod,size_t sample_size, unsigned int topN, 
                    unsigned int K, unsigned int K_UMI){
    int err = NO_ERROR;
    unsigned int idx; 
    double epsilon = 1e-40;

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

    for (unsigned int b = 0; b < K_UMI; b++){
        for (unsigned int k = 0; k < K; k++){
            idx = b * K + k;
            gamma[idx] /=  (eta[b]+epsilon);  // avoid /0
        }
        eta[b] /= sample_size;
    }

    return err;
}/* M_step */

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

    //[TODO] use the maginal likelihood to assign reads and barcodes                    
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