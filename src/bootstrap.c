#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void bootSample(size_t* bootCount, double* probs, size_t total, size_t countSize, size_t seed)
{

  int i, j;

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, seed);

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
