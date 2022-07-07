#ifndef __FFT_H__
#define __FFT_H__

#ifdef FFTW_LONG_DOUBLE
#	define FFTW_COMPLEX		fftwl_complex
#	define FFTW_PLAN		fftwl_plan
#	define FFTW_MALLOC		fftwl_malloc
#	define FFTW_PLAN_DFT_1D		fftwl_plan_dft_1d
#	define FFTW_PLAN_DFT_2D		fftwl_plan_dft_2d
#	define FFTW_EXECUTE		fftwl_execute
#	define FFTW_DESTROY_PLAN	fftwl_destroy_plan
#	define FFTW_FREE		fftwl_free
#elif defined FFTW_FLOAT
#	define FFTW_COMPLEX		fftwf_complex
#	define FFTW_PLAN		fftwf_plan
#	define FFTW_MALLOC		fftwf_malloc
#	define FFTW_PLAN_DFT_1D		fftwf_plan_dft_1d
#	define FFTW_PLAN_DFT_2D		fftwf_plan_dft_2d
#	define FFTW_EXECUTE		fftwf_execute
#	define FFTW_DESTROY_PLAN	fftwf_destroy_plan
#	define FFTW_FREE		fftwf_free
#else
#	define FFTW_COMPLEX		fftw_complex
#	define FFTW_PLAN		fftw_plan
#	define FFTW_MALLOC		fftw_malloc
#	define FFTW_PLAN_DFT_1D		fftw_plan_dft_1d
#	define FFTW_PLAN_DFT_2D		fftw_plan_dft_2d
#	define FFTW_EXECUTE		fftw_execute
#	define FFTW_DESTROY_PLAN	fftw_destroy_plan
#	define FFTW_FREE		fftw_free
#endif


#include <stdio.h>
#include <math.h>
#include <complex.h>	/* triggers use of C99 complex numbers */
#include <fftw3.h>	/* very important to use version 3 */

int fft(int n, FFTW_COMPLEX *f, int sign);
int fft2(int m, int n, FFTW_COMPLEX **f, int sign);

#endif
