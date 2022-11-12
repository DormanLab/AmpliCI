/**
 * @file fft.c
 * @author R. Maitra
 * @date 2005/08/26
 *
 * Calculates the FFT of complex periodic 1-, 2- and 3-dimensional arrays,
 * using FFTW, described at www.fftw.org.  Note that the FFT's are scaled
 * which means that applying the FFT and then its inverse provides the
 * sequence back.  This is a variation on the output of FFTW, which does
 * not scale.
 *
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"


/**
 * Calculate the forward and backward Fast Fourier Transform of a complex,
 * 1-dimensional sequence of period \p n.  It uses the weighting 1/n for
 * the forward transform and 1 for the inverse transform.
 *
 * @param n	periodicity (and size) of sequence
 * @param f	pointer to input sequence (overwritten)
 * @param sign	choose direction (FFTW_FORWARD|FFTW_BACKWARD)
 * @return	error status
 */
int fft(int n, FFTW_COMPLEX *f, int sign)
{
	FFTW_COMPLEX *in, *out;
	FFTW_PLAN p;
	int i;

	if (!abs(sign)) {
		fprintf(stderr, "ERROR: Please specify one of FFTW_FORWARD or "
			"FFTW_BACKWARD for forward/inverse transform\n");
		exit(1);
	}

	in = FFTW_MALLOC(n * sizeof *in);
	out = FFTW_MALLOC(n * sizeof *out);

	p = FFTW_PLAN_DFT_1D(n, in, out, sign, FFTW_MEASURE);
	for (i = 0; i < n; ++i)
		in[i] = f[i];

	FFTW_EXECUTE(p); /* repeat as needed */
	for (i = 0; i < n; ++i)
		if (sign == FFTW_FORWARD)
			f[i] = out[i] / n;
		else
			f[i] = out[i];

	FFTW_DESTROY_PLAN(p);
	FFTW_FREE(in);
	FFTW_FREE(out);
	return 0;
}

/**
 * Calculate foward and inverse Fast Fourier Transform of a complex
 * 2-dimensional m x n array.  See fft()
 *
 * @param m	first dimension
 * @param n	second dimension
 * @param f	input matrix
 * @param sign	forward or inverse
 * @return	error status
 */
int fft2(int m, int n, FFTW_COMPLEX **f, int sign)
{
	FFTW_COMPLEX *in, *out;
	FFTW_PLAN p;
	int i, j;

	if (!abs(sign)) {
		fprintf(stderr, "ERROR: Please specify one of FFTW_FORWARD or "
			"FFTW_BACKWARD for forward/inverse transform\n");
		exit(1);
	}

	in = FFTW_MALLOC(m * n * sizeof *in);
	out= FFTW_MALLOC(m * n * sizeof *out);

	p = FFTW_PLAN_DFT_2D(m, n, in, out, sign, FFTW_MEASURE);

	for (i = 0; i < m; ++i)
		for (j = 0; j < n; ++j)
			in[i * n + j] = f[i][j];

	FFTW_EXECUTE(p); /* repeat as needed */
	for (i = 0; i < m ; ++i)
		for (j = 0; j < n; ++j)
			if (FFTW_FORWARD)
				f[i][j] = out[i*n + j]/ m * n;
			else
				f[i][j] = out[i*n + j];
	FFTW_DESTROY_PLAN(p);
	FFTW_FREE(in);
	FFTW_FREE(out);
	return 0;
}
