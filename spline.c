#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void bootSample(size_t* bootCount, double* probs, size_t total, size_t countSize)
{

  int i, j;

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  size_t left;
  size_t runningTotal;
  size_t current;
  double correct;

  for (i = 0; i < 100; i++)
    {
      left = total;
      correct = 1;
      runningTotal = 0;
      current = 0;
      for (j = 0; j < countSize; j++)
	{
	  current = gsl_ran_binomial (r, probs[j] / correct, left);
	  runningTotal += current;
	  bootCount[countSize * i + j] = runningTotal;

	  left -= current;
	  correct -= probs[j];
	}
    }
  gsl_rng_free (r);
}



void splineFit(double* lambda, double* pi0, double* pi0Est, size_t length, int ncoeff)
{
  int nbreak = ncoeff - 2;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  gsl_vector *c;
  gsl_vector *x, *y;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq;

  /* allocate a cubic bspline workspace (k = 4) */
  bw = gsl_bspline_alloc(4, nbreak);
  B = gsl_vector_alloc(ncoeff);

  x = gsl_vector_alloc(length);
  y = gsl_vector_alloc(length);
  X = gsl_matrix_alloc(length, ncoeff);
  c = gsl_vector_alloc(ncoeff);
  cov = gsl_matrix_alloc(ncoeff, ncoeff);
  mw = gsl_multifit_linear_alloc(length, ncoeff);

  /* this is the data to be fitted */
  for (i = 0; i < length; ++i)
    {
      gsl_vector_set(x, i, lambda[i]);
      gsl_vector_set(y, i, pi0[i]);
    }

  /* use uniform breakpoints on [0, 15] */
  gsl_bspline_knots_uniform(lambda[0], lambda[length - 1], bw);

  /* construct the fit matrix X */
  for (i = 0; i < length; ++i)
    {
      double xi = gsl_vector_get(x, i);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      /* fill in row i of X */
      for (j = 0; j < ncoeff; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
        }
    }

  /* do the fit */
  gsl_multifit_linear(X, y, c, cov, &chisq, mw);

  //  output the smoothed curve
  double yerr;
  for (i = 0; i < length; i++)
    {
      gsl_bspline_eval(lambda[i], B, bw);
      gsl_multifit_linear_est(B, c, cov, pi0Est++, &yerr);
    }

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

}
