/**
 * @file bp_fft.c
 * @author K. S. Dorman
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
#include "fft.h"

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

int main(void) 
{
	FFTW_COMPLEX *dat;
	double eff = 0.8;
	double sam = 1.;//1.0;//0.2;
	int ngen = 3;
	int i, n = (1 << ngen) + 1;	/* pow(10, ngen) guarantees perfect
					 * precision, but support is practically
					 * much smaller
					 */
	double mean = 0, sum = 0, s;

	if (n < 2)
		n = 2;


	dat = FFTW_MALLOC(n * sizeof *dat);
	for (i = 0; i < n; ++i) {
		s = 2 * M_PI * i / n;
		dat[i] = pgf(ngen, cos(s) + sin(s) * I, eff, sam);
		//fprintf(stderr, "%f %fi\n", creal(dat[i]), cimag(dat[i]));
	}

	fft(n, dat, FFTW_FORWARD);

	for (i = 0; i < n; ++i) {
		mean += i * creal(dat[i]);
		sum += creal(dat[i]);
		printf("%d %f %fi\n", i, creal(dat[i]), cimag(dat[i]));
	}
	double mean1 = 1 + eff;
	printf("Mean: %f (expected mean: %f)\n", mean, pow(mean1, ngen) * sam);
  
	FFTW_FREE(dat);
	return 0;
}
