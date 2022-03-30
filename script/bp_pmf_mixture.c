/**
 * @file bp_pmf.c
 * @author Xiyu Peng, K. S. Dorman
 *
 * Demonstrate use of fft for branching process pmf.  Requires fft.c and fft.h
 * from fft directory.
 *
 * Compile as: varying precisions from double, long double, to float
gcc -Wall -pedantic -o bp_pmf bp_pmf.c fft.c -lfftw3 -lm
gcc -Wall -pedantic -o bp_pmf bp_pmf.c fft.c -DFFTW_LONG_DOUBLE -lfftw3l -lm
gcc -Wall -pedantic -o bp_pmf bp_pmf.c fft.c -DFFTW_FLOAT -lfftw3f -lm
 *
 * Learning objectives:
 * - use of FFTW library
 * - use of FFT to estimate PGFs
 * - the low precision of pow() and even powl()
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

//#define MATHLIB_STANDALONE 1
//#include <Rmath.h>
// #include <R.h>

#include "fft.h"

typedef struct {
	double epsilon;  // PCR error rate 
	double delta;   // sequence error rate
	unsigned int ncycles;
	double E;     //Efficency
	int input;   // user input these parameters 
} itr_para ;

/**
 * Probability generating function.
 *
 * I spent hours trying to figure out a problem with my initial implementation:
 * pow() even powl() is EXTREMELY low precision, which is why my solution here
 * uses explicit multiplication.
 *
 * @param ngen	number of generations
 * @param s	complex argument
 * @param eff	probability of replication
 * @param sam	sampling probability
 * @return	complex value of function
 */
FFTW_COMPLEX pgf(int ngen, FFTW_COMPLEX s, double eff, double sam)
{
	s = (1 - sam) + sam * s;
	for (int i = 0; i < ngen; ++i)
		s = (1 - eff) * s + eff * s * s;
	return s;
} /* pgf */


int fread_uints(FILE *fp, unsigned int *id, size_t n)
{
	id[0] = 0;
	for (unsigned int i = 1; i < n; ++i){  // n = high_thres + 1
		if (fscanf(fp, "%u", &id[i]) != 1) 
			id[i] = 0;
		//printf("%d %d \n", i, id[i]);
	}

	return 0;  // error status
} /* fread_uints */

/* P is the data and Q is the approximate */
double KL_div(double* P, double *Q, unsigned int n){
	double div = 0.;

	for(unsigned int i = 0; i < n; i++){
		if (Q[i] <= 0 ) Q[i] = 1e-50; // avoid numerical error
		if(P[i]> 0)
			div += P[i]* log(P[i]) - P[i] * log(Q[i]);
	}

	return div;
}

double KS_test_stat(double* P, double *Q, unsigned int n){
	double D = 0.;
	double cumP = 0.;
	double cumQ = 0.;
	double diff;

	for(unsigned int i = 0; i < n; i++){
			cumP+= P[i];
			cumQ+= Q[i];
			if (cumP > cumQ)
				diff = cumP- cumQ;
			else 
				diff = cumQ- cumP;
			if (diff > D)
				D = diff;	
	}

	return D;
}

double multi_ll(unsigned int *count, double *prob, unsigned int n){
	double ll = 0.;

	for(unsigned int i = 0; i < n; i++){
		if(prob[i]> 0)
			ll += count[i] * log(prob[i]);
	}

	return -ll;
}


double bp_pmf(double eff, unsigned int ngen, double sam, int*len, double**pmf){

	FFTW_COMPLEX *dat;


	int i, n = pow(2, ngen)+1;	/* pow(10, ngen) guarantees perfect
					 * precision, but support is practically
					 * much smaller
					 */
	*len = n+1;

	double mean = 0, sum = 0, s;

	if (n < 2)
		n = 2;

	/* allocate space for pmf */
	double *pm = calloc(n, sizeof (double));
	if(!pm)
		return 0;
	
	dat = FFTW_MALLOC(n * sizeof *dat);
	for (i = 0; i < n; ++i) {
		s = 2 * M_PI * i / n;
		dat[i] = pgf(ngen, cos(s) + sin(s) * I, eff, sam);
	}

	fft(n, dat, FFTW_FORWARD);

	for (i = 0; i < n; ++i) {
		mean += i * creal(dat[i]);
		sum += creal(dat[i]);
		pm[i] = creal(dat[i]);
		// printf("%d %f %fi\n", i, creal(dat[i]), cimag(dat[i]));
	}
	/* correct the bug */
	//pm[n] = pm[0];
	//pm[0] = 0.;
	//mean+= n * pm[n];
	
	double mean1 = 1 + eff;
	//printf("Mean: %f (expected mean: %f)\n", mean, pow(mean1, ngen) * sam);
  
	*pmf = pm;
	FFTW_FREE(dat);
	return 0;
}

/* branching process model with errors modeled */
double bp_err(unsigned int *obser_abun, unsigned int n, double E, 
		unsigned int ncycles,double epsilon, double delta,int output){

	int *ns;   // length of non-zero entries
	double **pmfs; // pmf of error free molecules

	double *R;  // P(R = i) pmf
	double *pmf_Y; // pmf of error reads 
	double *pmf_mix; // pmf of the mixture model
	// double *pmf_Z; // pmf of error reads with sequence error
	double *sens;
	double *prec;

	unsigned int Nplus1 = ncycles + 1;
	unsigned int Nplus2 = ncycles + 2;
	double Eplus1 = E + 1;

	ns = calloc(Nplus1, sizeof (int));
	if(!ns)
		return 0;

	pmfs = malloc(Nplus1 * sizeof (double*)); // N + 1 
	if(!pmfs)
		return 0;

	R = calloc(Nplus2, sizeof (double));
	if(!R)
		return 0;

	/* abundance distribution of error free reads */
	double Et = E * (1.0-epsilon);
	/* assume %1 reads would be read as an error. thus 0.99 sample rates */
	// sample rate could be 1-\epsilon

	/* cycle 0 */
	ns[0] = 2;
	double pfm0[2];
	pmfs[0] = &pfm0[0];
	pmfs[0][0] = 0.;
	pmfs[0][1] = 1.;

	/* cycle 1 to cycle N, error free pmf */

	for(unsigned int i = 1; i < ncycles; i++)
		bp_pmf(Et, i, 1., &ns[i], &pmfs[i]);
	bp_pmf(Et, ncycles, 1., &ns[ncycles], &pmfs[ncycles]);
	 

	pmf_Y = calloc(ns[ncycles], sizeof (double));
	if(!pmf_Y)
		return 0;

	sens = calloc(ns[ncycles], sizeof (double));
	if(!sens)
		return 0;

	prec = calloc(ns[ncycles], sizeof (double));
	if(!prec)
		return 0;

	pmf_mix = calloc(ns[ncycles], sizeof (double));
	if(!pmf_mix)
		return 0;
	
	/* error pmf */
	/* P(R = i) */
	double sumR = 0.;
	for(unsigned int i = 1; i < Nplus1; i++){
		R[i] = pow(Eplus1, i-1) * epsilon;  //PCR error 
		sumR += R[i];
	}
	R[Nplus1] = pow(Eplus1, ncycles) * delta;  //sequencing error
	sumR += R[Nplus1];
	
	//for(unsigned int i = ncycles; i > 0; i--)
	//	R[i] = R[i] - R[i-1];
	R[0] = 0;

	// check the proportions of errors in each cycle
	double sumR2 = 0.;
	for(unsigned int i = 0; i < Nplus2; i++){
		R[i] = R[i]/sumR;
		sumR2 += R[i];
		// printf("R: %f \n", R[i]);
	}

	// printf("sum_R: %f \n", sumR2);


	/* P(Y = y) */
	unsigned int lenN = ns[ncycles];
	for (unsigned int y = 0; y < lenN; ++y)
		for(unsigned int i = 1; i < Nplus1; i++)
			if((int) y < ns[ncycles-i])
				pmf_Y[y] += pmfs[ncycles-i][y] * R[i];
	pmf_Y[1] += R[Nplus1];

	//for(unsigned int i = 1; i < Nplus1; i++)
	//	fprintf(stderr, "pdf at 0: %f\n", pmfs[ncycles-i][0]);


	/* match the observe data */
	double emp_pmf[n];

	double sum_obser = 0;
	for (unsigned int y = 0; y < n; ++y)
		sum_obser += obser_abun[y];

	for (unsigned int y = 0; y < n; ++y)
		emp_pmf[y] = (double) obser_abun[y] / sum_obser;


	double cdf_Y = 0.;
	double cdf_X = 0.;
	for (unsigned int y = 0; y < lenN; ++y){
		//fprintf(stderr, "%f %f \n", pmfs[ncycles][y],pmf_Y[y]);
		cdf_Y += pmf_Y[y];  // P(X <= x | error )
		cdf_X += pmfs[ncycles][y];   // P(X <= x | true )
		// mean_Y += y* pmf_Y[y];
	}

	//fprintf(stderr, "cdfY_sum: %f\n", cdf_Y);
	//fprintf(stderr, "cdfX_sum: %f\n", cdf_X);


	/* consider the mixture model and measure the KL divergence */
	double lambdae = sumR + 1;
	lambdae = 1. - 1./lambdae ;
	fprintf(stderr, "lambdae: %f \n", lambdae);

	
	for (unsigned int y = 0; y < lenN; ++y)
		pmf_mix[y] = lambdae * pmf_Y[y] + (1- lambdae) * pmfs[ncycles][y]; 
	

	/* truncate the distribution */
	unsigned int hthres = (n > lenN)? lenN: n;
	unsigned int lthres = 1;
	unsigned int total_num = 0;

	double trun_cdf_mix = 0.;
	double trun_cdf_emp = 0.;
	for (unsigned int y = lthres; y < hthres; ++y){
		trun_cdf_mix += pmf_mix[y];
		trun_cdf_emp += emp_pmf[y];
		total_num += obser_abun[y];
	}

	sens[0] = 0;
	prec[0] = 0;
	for(unsigned int y = 1; y < hthres; ++y){
		sens[y] = pmfs[ncycles][y] + sens[y-1]; // P(X < = x)
		prec[y] = pmf_Y[y]+ prec[y-1]; // P(Y < = y)
	}
	
	
	for (unsigned int y = 0; y < hthres; ++y){
		pmf_mix[y] /= trun_cdf_mix;   // truncated pmf in [1,n]
		emp_pmf[y] /= trun_cdf_emp;
		sens[y] = sens[y] / sens[hthres-1];
		prec[y] = prec[y]/prec[hthres-1];
		if(output){ 
			if(!y){
				fprintf(stdout, "x, P(X = x | Z = 1), P(X = x | Z = 0), P_trun(X = x), P_emp(X = x), count, P(X > x | Z = 1), " 
								"P(X <= x | Z = 1), P(X <= x | Z = 0), P(X > x | Z = 0), "
								"P(X > x | Z = 1)+ P(X <= x | Z = 0) \n");
			}else{
				fprintf(stdout, "%i, %f, %f, %f, %f, %i, %f, %f, %f, %f, %f\n", y, pmfs[ncycles][y], pmf_Y[y], pmf_mix[y], 
					emp_pmf[y], obser_abun[y], 1-sens[y],sens[y], prec[y], 1-prec[y], 1-sens[y]+prec[y]);
			}
		}
	}

	double D = KS_test_stat(&emp_pmf[lthres], &pmf_mix[lthres], hthres-lthres);
	double ll = multi_ll(&obser_abun[lthres], &pmf_mix[lthres], hthres-lthres);  // -log likelihood
	double bic = 4 * log(total_num) + 2*ll;

	// fprintf(stderr, "D: %f; ll: %f\n", D, ll);
	
	free(ns);
	free(R);
	free(pmf_Y);
	for(unsigned int i = 1; i < Nplus1; i++)
		free(pmfs[i]);
	free(pmfs);
	free(pmf_mix);
	free(sens);
	free(prec);
	
	return D;
}

int PCR_para(unsigned int *obser_abun, unsigned int n, double *E_est, unsigned int *ncycles_est, 
				double epsilon0, double delta0, double *dist){

	double E_vec[20];
	unsigned int ncycles_vec[20];
	double dist_vec[400];

	double E_opt= 0., min_dist = INFINITY;
	unsigned int ncycles_opt = 1;

	
	E_vec[0] = 0.05;
	ncycles_vec[0] = 6;
	for(int i = 1; i < 20; ++i){
		E_vec[i] = E_vec[i-1]+0.05;
		ncycles_vec[i] = ncycles_vec[i-1] + 1;
	}

	unsigned int cnt = 0;
	for(int i = 0; i < 20; ++i){
		for (int j = 0; j < 16; ++j){
			double E = E_vec[i];
			double ncycles = ncycles_vec[j];   // error proportion 
			double dist = bp_err(obser_abun, (unsigned int) n, E, ncycles,epsilon0,delta0,0);
			dist_vec[cnt] = dist;
			fprintf(stderr, "[%f, %i], div: %f\n",  E_vec[i], ncycles_vec[j], dist_vec[cnt]);
			if(dist_vec[cnt] < min_dist){
				min_dist = dist_vec[cnt];
				E_opt =  E_vec[i];
				ncycles_opt = ncycles_vec[j];
			}
			cnt ++;
		}
	} 

	*E_est = E_opt;
	*ncycles_est = ncycles_opt;
	*dist = min_dist;

	return 0;
}

int err_para(unsigned int *obser_abun, unsigned int n, double E0, unsigned int ncycles0, 
				double *epsilon_est, double *delta_est, double *dist){
	
	double epsilon_vec[10];
	double delta_vec[10];
	double div_vec[100];

	double epsilon_opt= 0., min_div = INFINITY;
	double delta_opt = 0.;

	
	epsilon_vec[0] = 0.0005;
	delta_vec[0] = 0.005;
	for(int i = 1; i < 10; ++i){
		epsilon_vec[i] = epsilon_vec[i-1]+0.0005;
		delta_vec[i] = delta_vec[i-1] + 0.001;
	}

	unsigned int cnt = 0;
	for(int i = 0; i < 10; ++i){
		for (int j = 0; j < 10; ++j){

			double div = bp_err(obser_abun, (unsigned int) n, E0, ncycles0,epsilon_vec[i],delta_vec[j],0);
			div_vec[cnt] = div;
			fprintf(stderr, "[%f, %f], div: %f\n", epsilon_vec[i], delta_vec[j], div_vec[cnt]);
			if(div_vec[cnt] < min_div){
				min_div = div_vec[cnt];
				epsilon_opt =  epsilon_vec[i];
				delta_opt = delta_vec[j];
			}
			cnt ++;
		}
	} 

	*epsilon_est = epsilon_opt;
	*delta_est = delta_opt;
	*dist = min_div;

	return 0;
}


int para_est(unsigned int *obser_abun, unsigned int n, 
	double *E_est, unsigned int *ncycles_est, double *epsilon_est, double *delta_est){

	double E = 0.5;
	unsigned int ncycles = 15;
	double epsilon = 0.005;
	double delta = 0.01;
	double dist = 1;


	PCR_para(obser_abun, n, &E, &ncycles, epsilon, delta, &dist);
	fprintf(stderr, "FINAL: [%f, %i], div: %f\n", E, ncycles, dist);

	double dist_pre = dist + 100000;
	while( (dist_pre - dist)/dist > 1e-4 ){
		
		dist_pre = dist;
		err_para(obser_abun, n, E, ncycles, 
					&epsilon, &delta, &dist);
		fprintf(stderr, "FINAL: [%f, %f], div: %f\n", epsilon, delta, dist);

		if((dist_pre - dist)/dist < 1e-4)
			break;

		dist_pre = dist;
		PCR_para(obser_abun, n, &E, &ncycles, epsilon, delta, &dist);
		fprintf(stderr, "FINAL: [%f, %i], div: %f\n", E, ncycles, dist);
		
	}

	*E_est = E;
	*ncycles_est = ncycles;
	*epsilon_est = epsilon;
	*delta_est = delta;

	return 0;
}


int parse_options_simple(int argc, const char **argv, itr_para* para, const char** file_name, int *hthres){	
	int i, j;
	size_t n;
	int err = 0;
	char a;
	unsigned int hthre;

	for (i = 1; i < argc; i++) {

		n = strlen(argv[i]);
		if (n < 1)
			goto CMDLINE_ERROR;
		j = 0;
		while ((a = argv[i][++j]) == '-' && j < (int) n);
		switch(a) {
			case 't':
				if (i == argc - 1) {
					err = 1;
					goto CMDLINE_ERROR;
				} else {
					if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57){
							hthre = strtoul(argv[++i], NULL, 0);
							// fprintf(stderr, "%i\n",hthre);
							*hthres = (int) hthre;
						}
					break;	
				}
				break;
			case 'n':
				if (i == argc - 1) {
					err = 1;
					goto CMDLINE_ERROR;
				} else {
					if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57){
							para->ncycles = strtoul(argv[++i], NULL, 0);
							para->input = 1;
						}
					break;	
				}
				break;
			case 'd':
				if (i == argc - 1) {
					err = 1;
					goto CMDLINE_ERROR;
				} 
				para->delta = strtod(argv[++i], NULL);
				para->input = 1;
				break;
			case 'f':
				if (i == argc - 1){
					err = 1;
					goto CMDLINE_ERROR;
				}else{
				*file_name = argv[++i];
				}
				break;
			case 'e': 
				if (i == argc - 1) {
					err = 1;
					goto CMDLINE_ERROR;
				}
				if (!strncmp(&argv[i][j], "eff", 3)) {
					para->E = strtod(argv[++i], NULL);
					para->input = 1;
					break;
				}else{
					para->epsilon = strtod(argv[++i], NULL);
					para->input = 1;
					break;
				}
				break;	
			default:
				goto CMDLINE_ERROR;
		}
	}
	return err;

CMDLINE_ERROR:
	if(err > 0)
		fprintf(stderr, "INVALID_CMD_ARGUMENT !\n");
	return err;

} 


int main(int argc, const char **argv) {

	/* The file records the UMI observed abundance distribution. 
		Check sim8.2_raw_abun.txt as an example in the test folder */

	const char *file_name = NULL; //"sim8.2_raw_abun.txt";   //  

	/* truncation choose on the right tail */
	int hthres = 100; 

	/* set up */
    itr_para para = {
		.epsilon = 0.0045,
        .delta = 0.013,
		.ncycles = 10,
        .E = 0.5,
		.input = 0,
	};

	/* read from the command line */
	if(parse_options_simple(argc, argv, &para, &file_name, &hthres))
		return 0;

	/*  major parameter */
	double epsilon = para.epsilon;  // PCR error rate 
	double delta = para.delta;   // sequence error rate
	unsigned int ncycles = para.ncycles;
	double E = para.E;
	
	// hthres = 100; 
	fprintf(stderr, "INIT: [E: %f, ncycles: %i, epsilon: %f, delta: %f, H thres: %i]\n", E, ncycles, epsilon, delta, hthres);

	int n = hthres + 1;  // include 0
	unsigned int * obser_abun = NULL;
	obser_abun = calloc(n, sizeof (unsigned int));
	if(!obser_abun)
		return 0;

	/* read observed abundance <= hthres */
	FILE *fp = fopen(file_name, "r");
	if (!fp)
		return 0;
	fread_uints(fp, obser_abun,n);		
	fclose(fp);
	fprintf(stderr,"Finish reading %s \n", file_name);

	/* estimate the parameters */
	if(!para.input)
		para_est(obser_abun, (unsigned int) n, &E, &ncycles, &epsilon, &delta);

	/* run again to get the final distribution */
	double dist = bp_err(obser_abun, (unsigned int) n, E, ncycles, epsilon, delta, 1);

	fprintf(stdout, "FINAL: [E: %f, ncycles: %i, epsilon: %f, delta: %f], div: %f\n", E, ncycles, epsilon, delta, dist);
	fprintf(stderr, "FINAL: [E: %f, ncycles: %i, epsilon: %f, delta: %f], div: %f\n", E, ncycles, epsilon, delta, dist);

	free(obser_abun);
	return 0;
}
