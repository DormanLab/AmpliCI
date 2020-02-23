#ifndef __H_UTIL__
#define __H_UTIL__

#include <stdio.h>
#include <stddef.h>

int fread_size_ts(FILE *fp, size_t *idx, size_t n);
int fread_uints(FILE *fp, unsigned int *id, size_t n);
int read_uints(char const * const filename, unsigned int *id, size_t nreads);

#endif
