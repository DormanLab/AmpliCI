#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__

#include "error.h"

/**
 * Nelder-Mead errors.
 */
enum {
	OUT_OF_BAND = NUM_ERRORS + 1	/* outside band in banded alignment */
};

unsigned char **nwalign(unsigned char const * const s1, unsigned char const * const s2, 
	size_t len1, size_t len2, int score[4][4], int gap_p, int band, int ends_free, 
	double const *perr, int *err, size_t *alen);
int ana_alignment(unsigned char**aln, size_t alen, unsigned int rlen, unsigned int* nindels, 
					unsigned int *nmismatch, int ends_free,int dbg);

#endif
