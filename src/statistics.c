#include "statistics.h"
#include "error.h"
double aic(double ll, size_t k) {
	return (2*k - 2*ll);
} /* aic */

double bic(double ll, size_t k, size_t n) {
	return (log(n)*k - 2*ll);
} /* bic */

/**
 * Compute the pmf for random variable X ~ Poisson binomial(n, perr).
 *
 * modified from code in STAT580
 *
 * @param n	number of trials, length of perr
 * @parma k	value of X
 * @param perr	probabilities of success
 * @return	probability
 */
double dpoisbind(unsigned int n, unsigned int k, double *perr)
{
	double pnj[k + 1];

	/* initialize: p1k */
	pnj[0] = 1 - perr[0];
	if (k)
		pnj[1] = perr[0];

	for (unsigned int j = 2; j <= n; ++j) {
		if (k) {	/* pjk */
			if (k >= j)
				pnj[j] = perr[j - 1] * pnj[j - 1];
			for (unsigned int i = k >= j ? j - 1 : k; i > 0; --i)
				pnj[i] = perr[j - 1] * pnj[i - 1]
						+ (1 - perr[j - 1]) * pnj[i];
		}
		pnj[0] = (1 - perr[j-1]) * pnj[0];
//		fprintf(stderr, "%e\n", pnj[k]);
	}
	
	return pnj[k];
}/* dpoisbind */


/**
 * Probability distribution of X ~ Poisson binomial(n, perr).  Compute
 * Pr(X <= k) or Pr(X > k) as per \par upper_tail.
 *
 * [TODO] need more efficient algorithm
 *
 * @param k		compute Pr(X <= k)
 * @param n		number of trials, length of perr
 * @param perr		probabilities of success
 * @param upper_tail	compute Pr(X>k)
 * @return		desired probability
 */
double ppoisbin(int k, unsigned int n, double *perr, int upper_tail)
{
	double prob = 0.;
	
	if (k < 0) {
		//mmessage(INFO_MSG, NO_ERROR, "1\n");
		if(upper_tail)
			return 1.;
		else
			return 0.;
	} else if ((unsigned int) k >= n) {
		//mmessage(INFO_MSG, NO_ERROR, "0\n");
		if (upper_tail) 
			return prob;
		else
			return 1 - prob;
	}

	/* [KSD] You should compute left tail probability if k is small. */
	/* currently most k are very large, close to n */
	for (unsigned int i = n; i > (unsigned int) k; i--)
		prob += dpoisbind(n, i, perr);

	if (upper_tail) 
		return prob;
	else
		return 1 - prob;

}/* ppoisbin */
