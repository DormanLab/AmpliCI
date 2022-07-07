#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "align.h"
#include "error.h"

/**
 * Needleman-Wunsch alignment.
 *
 * @param s1	first sequence
 * @param s2	second sequence
 * @param len1	length of first sequence
 * @param len2	length of first sequence
 * @param score	scores
 * @param gap_p
 * @param band	band, -1 for no band
 * @param ends_free	ends-free alignment
 * @param perr	error probability in read (second sequence)
 * @param alen	pointer to alignment length
 * @return	alignment
 */
unsigned char **nwalign(unsigned char const * const s1, unsigned char const * const s2,
	size_t len1, size_t len2, int score[4][4], int gap_p, int band,
	int ends_free, double const *perr, int *err, size_t *alen) {
	static size_t nnw = 0;
	size_t i, j;
	int l, r;   // BUG here
	size_t iband = band >= 0 ? band : 0;
	double diag, left, up;

	*err = NO_ERROR;

	unsigned int nrow = len1 + 1;
	unsigned int ncol = len2 + 1;

	//int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
	double *d = (double *) malloc(nrow * ncol * sizeof(double)); //E
	int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
	if (d == NULL || p == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "d & p");
		return NULL;
	}

	// Fill out left columns of d, p.
	for (i = 0; i <= len1; i++) {
		d[i*ncol] = ends_free ? 0 : i * gap_p; // ends-free gap
		p[i*ncol] = 3;
	}

	// Fill out top rows of d, p.
	for (j = 0; j <= len2; j++) {
		d[j] = ends_free ? 0 : j * gap_p; // ends-free gap
		p[j] = 2;
	}

	// Calculate left/right-bands in case of different lengths
	size_t lband, rband;
	if (len2 > len1) {
		lband = iband;
		rband = iband + len2 - len1;
	} else if (len1 > len2) {
		lband = iband + len1 - len2;
		rband = iband;
	} else {
		lband = iband;
		rband = iband;
	}

	// Fill out band boundaries of d.
	if (band >= 0 && (iband < len1 || iband < len2)) {
		for (i = 0; i <= len1; i++) {
			if ((int) i - (int) lband - 1 >= 0)
				d[i*ncol + i - lband - 1] = -9999;
			if (i + rband + 1 <= len2)
				d[i*ncol + i + rband + 1] = -9999;
		}
	}

	// Fill out the body of the DP matrix.
	for (i = 1; i <= len1; i++) {
		if (band >= 0) {
			l = i - lband;
			if (l < 1)
				l = 1;
			r = i + rband;
			if (r > (int) len2)
				r = len2;
		} else {
			l = 1;
			r = len2;
		}

		for (j = l; (int) j <= r; j++) {
			// Score for the left move.
			if (i == len1)
				left = d[i*ncol + j - 1]
					+ (ends_free ? 0 : gap_p); // Ends-free gap.
			else
				left = d[i*ncol + j - 1] + gap_p;

			// Score for the up move.
			if (j == len2)
				up = d[(i-1)*ncol + j]
					+ (ends_free ? 0 : gap_p); // Ends-free gap.
			else
				up = d[(i-1)*ncol + j] + gap_p;

			// Score for the diagonal move.
			diag = d[(i-1)*ncol + j-1]
				+ (perr ? perr[j-1] : 1.) * score[(int) s1[i-1]][(int) s2[j-1]];

			// Break ties and fill in d,p.
			if (up >= diag && up >= left) {
				d[i*ncol + j] = up;
				p[i*ncol + j] = 3;
			} else if (left >= diag) {
				d[i*ncol + j] = left;
				p[i*ncol + j] = 2;
			} else {
				d[i*ncol + j] = diag;
				p[i*ncol + j] = 1;
			}
		}
	}

	unsigned char *al0 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	unsigned char *al1 = (unsigned char *) malloc((len1+len2) * sizeof(unsigned char));
	if (al0 == NULL || al1 == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al0 & al1");
		return NULL;
	}

	// Trace back over p to form the alignment.
	size_t len_al = 0;
	i = len1;
	j = len2;

	//for (int ii = 0; ii < len1; ii++)
	//	mmessage(INFO_MSG, NO_ERROR, "d=%f\n", d[ii*ncol + ii]);

	while ( i > 0 || j > 0 ) {
	//	mmessage(INFO_MSG, NO_ERROR, "(%i, %i): p=%i, d=%f\n", i, j, p[i*ncol + j], d[i*ncol + j]);
		switch ( p[i*ncol + j] ) {
			case 1:
				al0[len_al] = s1[--i];
				al1[len_al] = s2[--j];
				break;
			case 2:
				al0[len_al] = '-';
				al1[len_al] = s2[--j];
				break;
			case 3:
				al0[len_al] = s1[--i];
				al1[len_al] = '-';
				break;
			default:
				*err = OUT_OF_BAND;
				mmessage(WARNING_MSG, *err,
					"NW alignment out of range");
				return NULL;
		}
		len_al++;
	}
	
	// Allocate memory to alignment strings.
	unsigned char **al = (unsigned char **) malloc( 2 * sizeof(unsigned char *) ); //E
	if (al == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al");
		return NULL;
	}

	al[0] = (unsigned char *) malloc(len_al); //E
	al[1] = (unsigned char *) malloc(len_al); //E
	if (al[0] == NULL || al[1] == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "al[]");
		return NULL;
	}

	// Reverse the alignment strings (since traced backwards).
	for (i = 0 ; i < len_al ; i++) {
		al[0][i] = al0[len_al-i-1];
		al[1][i] = al1[len_al-i-1];
	}

	// Free allocated memory
	free(d);
	free(p);
	free(al0);
	free(al1);

	*alen = len_al;

	nnw++;
	return al;
} /* nwalign */

/**
 * count num of indels and mismatches of Needleman-Wunsch alignment.
 *
 * @param aln	   Needleman-Wunsch alignment result
 * @param alen	   Alignment length 
 * @param rlen 		Reads (s2) length 
 * @param nindels   pointer to Number of indels 
 * @param nmismatch	pointer to Number of mismatches 
 * @param dbg	debug information
 * @return	0
 */
int ana_alignment(unsigned char**aln, size_t alen, unsigned int rlen, unsigned int* nindels, 
					unsigned int *nmismatch, int ends_free, int dbg){

		int fxn_debug = dbg;

		unsigned int nmis=0, nins = 0,ndel=0, nind = 0;
		if (aln) {
			for (size_t j= 0 ;j < alen; j++) {

				unsigned int j1 = j - nins;   // pos idx of hap
				unsigned int j2 = j - ndel;   // pos idx of read

				/* gaps in the end */
				if (j2 >= rlen || j1 >= rlen) // gaps in the end
					break;

				if (aln[0][j] == '-') {
					nins++;
					if (j == 0)
						nind += ends_free ? 0: 1;
					else if (aln[0][j-1] != '-')
						nind++;
					continue;
				}
	
				if (aln[1][j] == '-') {
					ndel++;
					if (j == 0)
						nind += ends_free ? 0: 1;
					else if (aln[1][j-1] != '-')
						nind++;
					continue;
				}
	
				if (aln[1][j] != aln[0][j])
					nmis++;
			}
		}

		debug_msg(DEBUG_III, fxn_debug, "num of indels: %i; num of "
			"mismatch: %i\n", nind, nmis);

		*nmismatch = nmis;
		*nindels = nind;

		return NO_ERROR;

}
