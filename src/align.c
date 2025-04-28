/**
 * @file align.c
 *
 * Banded Needleman-Wunsch alignments.  Semi-global alignments are achieved by
 * setting ends_free=1, however all the alignments conducted by this code
 * assume 3' gap penalty is 0 because Illumina reads are truncated.
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "align.h"
#include "error.h"

/**
 * Needleman-Wunsch alignment.
 *
 * @param s1		first sequence
 * @param s2		second sequence
 * @param len1		length of first sequence
 * @param len2		length of first sequence
 * @param score		scores
 * @param gap_p		gap penalty
 * @param band		band, -1 for no band
 * @param ends_free	ends-free alignment
 * @param perr		error probability in read (second sequence)
 * @param alen		pointer to alignment length
 * @param asc		alignment score
 * @return		alignment
 */
unsigned char **nwalign(unsigned char const * const s1, unsigned char const * const s2,
	size_t len1, size_t len2, int score[4][4], int gap_p, int band,
	int ends_free, double const *perr, int *err, size_t *alen, double *asc)
{
	int fxn_debug = ABSOLUTE_SILENCE;//DEBUG_I;//DEBUG_II;//
	static size_t nnw = 0;
	size_t i, j;
	int l, r;   // BUG here
	size_t iband = band >= 0 ? band : 0;
	double diag, left, up;
	int ends_free_3prime = 1;

	debug_msg(DEBUG_I, fxn_debug, "gap=%i band=%i semi=%i scores=\n", gap_p,
							band, ends_free);
	debug_msg_cont(DEBUG_I, fxn_debug, "%2i %2i %2i %2i\n%2i %2i %2i %2i\n"
		"%2i %2i %2i %2i\n%2i %2i %2i %2i\n",
		score[0][0], score[0][1], score[0][2], score[0][3],
		score[1][0], score[1][1], score[1][2], score[1][3],
		score[2][0], score[2][1], score[2][2], score[2][3],
		score[3][0], score[3][1], score[3][2], score[3][3]);

	*err = NO_ERROR;

	unsigned int nrow = len1 + 1;
	unsigned int ncol = len2 + 1;

	//int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
	double *d = calloc(nrow * ncol, sizeof(double)); //E
	int *p = calloc(nrow * ncol, sizeof(int)); //E
	if (d == NULL || p == NULL) {
		*err = MEMORY_ALLOCATION;
		mmessage(ERROR_MSG, *err, "d & p");
		return NULL;
	}

	// Fill out left column of d, p.
	d[0] = 0;
	for (i = 1; i <= len1; i++) {
		d[i*ncol] += ends_free ? 0 : gap_p; // ends-free gap
		p[i*ncol] = 3;
	}
	debug_msg(DEBUG_II, fxn_debug, "first col: %f\n", d[ncol]);

	// Fill out top row of d, p.
	for (j = 1; j <= len2; j++) {
		d[j] += ends_free ? 0 : gap_p; // ends-free gap
		p[j] = 2;
	}
	debug_msg(DEBUG_II, fxn_debug, "first row: %f\n", d[1]);

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
					+ (ends_free_3prime ? 0 : gap_p); // Ends-free gap.
			else
				left = d[i*ncol + j - 1] + gap_p;

			// Score for the up move.
			if (j == len2)
				up = d[(i-1)*ncol + j]
					+ (ends_free_3prime ? 0 : gap_p); // Ends-free gap.
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
			debug_msg(DEBUG_II, fxn_debug, "d[%i][%i]: %f (%.0f %.0f %.0f)\n", i, j, d[i*ncol + j], left, diag, up);
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
		//mmessage(INFO_MSG, NO_ERROR, "(%i, %i): p=%i, d=%f\n", i, j, p[i*ncol + j], d[i*ncol + j]);
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

	if (asc)
		*asc = d[len1*ncol + len2];

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
 * Count num of indels and mismatches of Needleman-Wunsch alignment.
 *
 * @param aln		Needleman-Wunsch alignment result
 * @param alen		Alignment length
 * @param rlen		Reads (s2) length
 * @param nindels	pointer to number of indel events
 * @param cnt_indels	pointer to number of indels
 * @param cnt_5prime	pointer to number of 5' indels
 * @param nmismatch	pointer to Number of mismatches
 * @param dbg		debug information
 * @return		error status
 */
int ana_alignment(unsigned char **aln, size_t alen, unsigned int rlen,
	unsigned int *nindels, unsigned int *cnt_indels,
	unsigned int *cnt_5prime_indels, unsigned int *nmismatch,
				int ends_free, int dbg)
{

	int fxn_debug = dbg;//DEBUG_I;//DEBUG_III;//
	unsigned int nmis = 0, nins = 0, ndel = 0, nind = 0;
	unsigned int n_5prime_ins = 0, n_5prime_del = 0;
	unsigned char started[2] = {0,0};

	debug_msg(DEBUG_I, fxn_debug, "sg=%i\n", ends_free);

	if (aln) {
		for (size_t j= 0 ;j < alen; j++) {

			unsigned int j1 = j - nins;   // pos idx of hap
			unsigned int j2 = j - ndel;   // pos idx of read

			/* gaps in the end */
			if (j2 >= rlen || j1 >= rlen) // gaps in the 3' end
				break;

			if (aln[0][j] == '-') {
				nins++;
				started[1] = 1;
				if (!started[0])
					++n_5prime_ins;
				if (j == 0)
					nind += ends_free ? 0: 1;
				else if (aln[0][j-1] != '-')
					nind++;
				continue;
			}

			if (aln[1][j] == '-') {
				ndel++;
				started[0] = 1;
				if (!started[1])
					++n_5prime_del;
				if (j == 0)
					nind += ends_free ? 0: 1;
				else if (aln[1][j-1] != '-')
					nind++;
				continue;
			}

			started[0] = started[1] = 1;

			if (aln[1][j] != aln[0][j])
				nmis++;
		}
	}

	debug_msg(DEBUG_III, fxn_debug, "id=%i (total=%i 5prime=%i) mm=%i\n",
		nind, nins + ndel, n_5prime_ins + n_5prime_del, nmis);

	*nmismatch = nmis;
	*nindels = nind;
	if (cnt_indels)
		*cnt_indels = nins + ndel;
	if (cnt_5prime_indels)
		*cnt_5prime_indels = n_5prime_ins + n_5prime_del;

	return NO_ERROR;

} /* ana_alignment */
