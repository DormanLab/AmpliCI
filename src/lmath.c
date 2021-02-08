#include <math.h>
#include "lmath.h"

/**
 * Compute squared Euclidean distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return squared distance between two vectors
 */
double sq_euc_dis(double *x, double *y, size_t p)
{
	double sum = 0;
	size_t i;

        for (i = 0; i < p; i++)
                sum += SQ(x[i] - y[i]);
        return sum;
} /* sq_euc_dis */

/**
 * Compute Euclidian distance between two p-vectors.
 *
 * @param x vector of length p
 * @param y vector of length p
 * @param p length of vector
 * @return distance between two vectors
 */
double euc_dis(double *x, double *y, size_t p)
{
	return sqrt(sq_euc_dis(x, y, p));
} /* euc_dis */

/**
 * Compute Hamming distance between two p-vectors.
 *
 * @param x	character vector of length p
 * @param y	character vector of length p
 * @param p	length of vector
 * @return	Hamming distance between two vectors
 */
size_t hamming_char_dis(char *x, char *y, size_t p)
{
	size_t hd = 0;
	for (size_t i = 0; i < p; ++i)
		hd += x[i] != y[i];
	return hd;
} /* hamming_char_dis */

/**
 * Compute Hamming distance between two p-vectors.
 *
 * @param x	character vector of length p
 * @param y	character vector of length p
 * @param p	length of vector
 * @return	Hamming distance between two vectors
 */
unsigned int hamming_uchar_dis(unsigned char *x, unsigned char *y, unsigned int p)
{
	unsigned int hd = 0;
	for (unsigned int i = 0; i < p; ++i)
		hd += x[i] != y[i];
	return hd;
} /* hamming_char_dis */

/**
 * Determine if given element is in a vector.
 * 
 * @param array	pointer to vector of numbers
 * @param num	element to find in vector
 * @param l	length of vector
 * @return	1/0 for true/false
 **/
int find_uint(unsigned int *array, unsigned int num, size_t l)
{
	for (size_t i = 0 ; i < l; i++)
		if (num == array[i])
			return 1;
	return 0;
}/* find_uint */
