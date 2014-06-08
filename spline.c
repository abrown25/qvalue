#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

double spline_fit(double* lambda, double* pi0, size_t length)
{
  const size_t ncoeffs = 12;
  const size_t nbreak = 10;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  double dy;
  gsl_vector *c;
  gsl_vector *x, *y;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq;

  /* allocate a cubic bspline workspace (k = 4) */
  bw = gsl_bspline_alloc(4, nbreak);
  B = gsl_vector_alloc(ncoeffs);

  x = gsl_vector_alloc(length);
  y = gsl_vector_alloc(length);
  X = gsl_matrix_alloc(length, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(length, ncoeffs);

  /* this is the data to be fitted */
  for (i = 0; i < length; ++i)
    {
      gsl_vector_set(x, i, lambda[i]);
      gsl_vector_set(y, i, pi0[i]);
    }

  /* use uniform breakpoints on [0, 15] */
  gsl_bspline_knots_uniform(0.0, 0.9, bw);

  /* construct the fit matrix X */
  for (i = 0; i < length; ++i)
    {
      double xi = gsl_vector_get(x, i);

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);

      /* fill in row i of X */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(X, i, j, Bj);
        }
    }

  /* do the fit */
  gsl_multifit_linear(X, y, c, cov, &chisq, mw);

  //  output the smoothed curve
  double xi, yi, yerr;
  gsl_bspline_eval(0.9, B, bw);
  gsl_multifit_linear_est(B, c, cov, &yi, &yerr);

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

  return yi;
}
