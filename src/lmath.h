#ifndef __MATH_TEST_H__
#define __MATH_TEST_H__

#include <stddef.h>
#include <math.h>

/**
 * Type of vectorization, by row or column.
 */
enum {
	ROW_ORDER,
	COLUMN_ORDER
};

/**
 * Simple mathematical functions.
 */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SQ(x) ((x) * (x))

double sq_euc_dis(double *x, double *y, size_t n);
double euc_dis(double *x, double *y, size_t n);
size_t hamming_char_dis(char *x, char *y, size_t n);
int find_uint(unsigned int *array, unsigned int num, size_t l);

#endif
