#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <math.h>
#include <stddef.h>

double aic(double ll, size_t k);
double bic(double ll, size_t k, size_t n);

/* possion binomial distribution */
double dpoisbind(unsigned int n, unsigned int k, double perr[]);
double ppoisbin(int k, unsigned int n, double *perr, int upper_tail);



#endif
