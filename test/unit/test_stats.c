/* Unit tests for the random number generators and the distribution
   functions in stats.c, and for the aLRT statistic conversion in alrt.c. */

#include "tinytest.h"
#include "utilities.h"
#include "stats.h"
#include "alrt.h"

#define N_DRAWS 200000

static double mean_of(double (*draw)(void *), void *arg, int n, double *var)
{
  double s = 0, s2 = 0;
  for (int i = 0; i < n; i++)
  {
    double x = draw(arg);
    s += x;
    s2 += x * x;
  }
  if (var) *var = s2 / n - (s / n) * (s / n);
  return s / n;
}

/* ------------------------------------------------------------ uniform */

TEST(random, uni_in_unit_interval)
{
  for (int i = 0; i < N_DRAWS; i++)
  {
    phydbl u = Uni();
    CHECK(u >= 0.0 && u <= 1.0);
  }
}

static double draw_uni(void *unused) { (void)unused; return Uni(); }

TEST(random, uni_moments)
{
  double var, m = mean_of(draw_uni, NULL, N_DRAWS, &var);
  CHECK_NEAR(m, 0.5, 0.005);
  CHECK_NEAR(var, 1.0 / 12.0, 0.002);
}

TEST(random, rand_int_within_bounds_and_covers_them)
{
  int seen[11] = {0};
  for (int i = 0; i < N_DRAWS; i++)
  {
    int k = Rand_Int(-5, 5);
    CHECK(k >= -5 && k <= 5);
    seen[k + 5]++;
  }
  for (int k = 0; k < 11; k++) CHECK(seen[k] > 0);
}

TEST(random, rand_int_single_value_range)
{
  for (int i = 0; i < 100; i++) CHECK_EQ_INT(Rand_Int(7, 7), 7);
}

TEST(random, rand_int_is_uniform)
{
  /* Chi-square goodness of fit over 10 bins, 9 df: 99.9% quantile 27.88. */
  int    counts[10] = {0};
  double chi2       = 0, expect = N_DRAWS / 10.0;
  for (int i = 0; i < N_DRAWS; i++) counts[Rand_Int(0, 9)]++;
  for (int k = 0; k < 10; k++) chi2 += pow(counts[k] - expect, 2) / expect;
  CHECK(chi2 < 27.88);
}

/* Bootstrap draws Rand_Int(0, n_sites - 1). With RAND_MAX = 32767 (MinGW,
   the Windows build) sites beyond 32,767 can never be drawn. This passes on
   platforms where RAND_MAX is 2^31 - 1. */
TEST(random, rand_int_reaches_large_values)
{
  int max = 0;
  for (int i = 0; i < N_DRAWS; i++)
  {
    int k = Rand_Int(0, 99999);
    if (k > max) max = k;
  }
  CHECK(max > 32767);
}

TEST(random, permutate_returns_a_permutation)
{
  for (int rep = 0; rep < 100; rep++)
  {
    int  n = 1 + rep % 17;
    int *p = Permutate(n);
    int  seen[17] = {0};
    for (int i = 0; i < n; i++)
    {
      CHECK(p[i] >= 0 && p[i] < n);
      seen[p[i]]++;
    }
    for (int i = 0; i < n; i++) CHECK_EQ_INT(seen[i], 1);
    Free(p);
  }
}

/* Permutate swaps x[i] with x[Rand_Int(0, len-1)] instead of
   x[Rand_Int(i, len-1)] (Fisher-Yates), which biases the permutations: for
   n = 3 some orders have probability 4/27 and others 5/27. */
KNOWN_BUG(random, permutate_is_uniform,
          "stats.c:4659 Permutate is not a Fisher-Yates shuffle")
{
  const int n_rep = 60000;
  int       counts[6] = {0};
  double    chi2 = 0, expect = n_rep / 6.0;
  for (int rep = 0; rep < n_rep; rep++)
  {
    int *p = Permutate(3);
    counts[p[0] * 2 + (p[1] > p[2])]++;
    Free(p);
  }
  /* 5 df, 99.9% quantile 20.52 */
  for (int k = 0; k < 6; k++) chi2 += pow(counts[k] - expect, 2) / expect;
  CHECK(chi2 < 20.52);
}

/* ------------------------------------------------- other generators */

static double draw_gamma(void *a) { double *p = a; return Rgamma(p[0], p[1]); }

TEST(random, rgamma_moments)
{
  double params[][2] = {{0.3, 2.0}, {1.0, 1.0}, {2.5, 0.4}, {10.0, 0.1}};
  for (int i = 0; i < 4; i++)
  {
    double var, m = mean_of(draw_gamma, params[i], N_DRAWS, &var);
    double em = params[i][0] * params[i][1];
    double ev = params[i][0] * params[i][1] * params[i][1];
    CHECK_NEAR(m, em, 0.02 * em);
    CHECK_NEAR(var, ev, 0.06 * ev);
  }
}

TEST(random, rgamma_positive)
{
  for (int i = 0; i < N_DRAWS; i++) CHECK(Rgamma(0.1, 1.0) >= 0.0);
}

static double draw_exp(void *a) { return Rexp(*(double *)a); }

TEST(random, rexp_mean)
{
  double lambda = 2.5;
  CHECK_NEAR(mean_of(draw_exp, &lambda, N_DRAWS, NULL), 1.0 / lambda, 0.01);
}

static double draw_norm(void *a) { double *p = a; return Rnorm(p[0], p[1]); }

TEST(random, rnorm_moments)
{
  double p[2] = {3.0, 2.0}, var;
  double m = mean_of(draw_norm, p, N_DRAWS, &var);
  CHECK_NEAR(m, 3.0, 0.03);
  CHECK_NEAR(var, 4.0, 0.08);
}

/* Used as Bayesian-bootstrap site weights: they sum to n, not 1. */
TEST(random, dirichlet_weights_are_positive_and_sum_to_n)
{
  for (int rep = 0; rep < 50; rep++)
  {
    int     n = 10 + rep;
    phydbl *w = Dirichlet(n);
    double  s = 0;
    for (int i = 0; i < n; i++)
    {
      CHECK(w[i] > 0.0);
      s += w[i];
    }
    CHECK_NEAR(s, (double)n, 1e-8);
    free(w);
  }
}

/* ------------------------------------------------- special functions */

TEST(special, ln_gamma)
{
  CHECK_NEAR(LnGamma(1.0), 0.0, 1e-10);
  CHECK_NEAR(LnGamma(5.0), log(24.0), 1e-10);
  CHECK_NEAR(LnGamma(0.5), 0.5 * log(M_PI), 1e-10);
  CHECK_NEAR(LnGamma(10.5), lgamma(10.5), 1e-9);
}

TEST(special, ln_fact_and_binomial)
{
  CHECK_NEAR(LnFact(0), 0.0, 1e-12);
  CHECK_NEAR(LnFact(5), log(120.0), 1e-10);
  CHECK_NEAR(Bico(5, 2), 10.0, 1e-9);
  CHECK_NEAR(Bico(10, 0), 1.0, 1e-9);
  CHECK_NEAR(Bico(20, 10), 184756.0, 1e-4);
}

TEST(special, incomplete_gamma)
{
  /* Regularised lower incomplete gamma P(a, x). */
  CHECK_NEAR(IncompleteGamma(1.0, 1.0, LnGamma(1.0)), 1.0 - exp(-1.0), 1e-9);
  CHECK_NEAR(IncompleteGamma(2.0, 3.0, LnGamma(3.0)), 1.0 - 5.0 * exp(-2.0), 1e-9);
  CHECK_NEAR(IncompleteGamma(0.0, 2.0, LnGamma(2.0)), 0.0, 1e-12);
}

TEST(special, normal_distribution)
{
  CHECK_NEAR(Pnorm(0.0, 0.0, 1.0), 0.5, 1e-9);
  CHECK_NEAR(Pnorm(1.959964, 0.0, 1.0), 0.975, 1e-6);
  CHECK_NEAR(Dnorm(0.0, 0.0, 1.0), 1.0 / sqrt(2 * M_PI), 1e-12);
  CHECK_NEAR(PointNormal(0.975), 1.959964, 1e-5);
  CHECK_NEAR(PointNormal(0.5), 0.0, 1e-6);
}

TEST(special, chi_square_quantiles)
{
  CHECK_NEAR(PointChi2(0.95, 1), 3.841459, 1e-4);
  CHECK_NEAR(PointChi2(0.99, 2), 9.210340, 1e-4);
  CHECK_NEAR(PointChi2(0.5, 10), 9.341818, 1e-4);
}

TEST(special, gamma_density)
{
  /* Dgamma(x, shape, scale) */
  CHECK_NEAR(Dgamma(1.0, 1.0, 1.0), exp(-1.0), 1e-10);
  CHECK_NEAR(Dgamma(2.0, 3.0, 0.5), pow(2.0, 2) * exp(-4.0) / (2.0 * pow(0.5, 3)), 1e-10);
  CHECK_NEAR(Pgamma(1.0, 1.0, 1.0), 1.0 - exp(-1.0), 1e-8);
}

/* --------------------------------------------------- discrete gamma */

TEST(discrete_gamma, yang_1994_table)
{
  /* Yang (1994) J. Mol. Evol. 39:306, alpha = 0.5, 4 categories (means). */
  phydbl freq[4], rate[4];
  DiscreteGamma(freq, rate, 0.5, 0.5, 4, 0);
  CHECK_NEAR(rate[0], 0.0334, 1e-4);
  CHECK_NEAR(rate[1], 0.2519, 1e-4);
  CHECK_NEAR(rate[2], 0.8203, 1e-4);
  CHECK_NEAR(rate[3], 2.8944, 1e-4);
  for (int i = 0; i < 4; i++) CHECK_NEAR(freq[i], 0.25, 1e-12);
}

TEST(discrete_gamma, rates_have_mean_one_and_increase)
{
  double alphas[] = {0.05, 0.3, 1.0, 4.0, 50.0};
  int    ks[]     = {2, 4, 8, 16};
  for (int a = 0; a < 5; a++)
    for (int k = 0; k < 4; k++)
      for (int median = 0; median < 2; median++)
      {
        phydbl freq[16], rate[16];
        double mean = 0;
        DiscreteGamma(freq, rate, alphas[a], alphas[a], ks[k], median);
        for (int i = 0; i < ks[k]; i++)
        {
          mean += rate[i] / ks[k];
          if (i) CHECK(rate[i] > rate[i - 1]);
        }
        CHECK_NEAR(mean, 1.0, 1e-6);
      }
}

TEST(discrete_gamma, one_category_is_rate_one)
{
  phydbl freq[1], rate[1];
  DiscreteGamma(freq, rate, 0.7, 0.7, 1, 0);
  CHECK_NEAR(rate[0], 1.0, 1e-12);
  CHECK_NEAR(freq[0], 1.0, 1e-12);
}

TEST(discrete_gamma, large_alpha_tends_to_equal_rates)
{
  phydbl freq[4], rate[4];
  DiscreteGamma(freq, rate, 1e4, 1e4, 4, 0);
  for (int i = 0; i < 4; i++) CHECK_NEAR(rate[i], 1.0, 0.03);
}

/* --------------------------------------------------------- aLRT */

TEST(alrt, statistic_to_probability_range_and_monotonicity)
{
  phydbl prev = -1.0;
  for (double x = 0.0; x < 20.0; x += 0.001)
  {
    phydbl p = Statistics_To_Probabilities(x);
    CHECK(p > 0.0 && p <= 1.0);
    CHECK(p >= prev - 1e-12);
    prev = p;
  }
}

TEST(alrt, statistic_to_probability_table_points)
{
  /* Interpolated chi2 table v, then ((1+v)/2)^3. */
  CHECK_NEAR(Statistics_To_Probabilities(3.841), pow(0.975, 3), 1e-9);
  CHECK_NEAR(Statistics_To_Probabilities(6.635), pow(0.995, 3), 1e-9);
  CHECK_NEAR(Statistics_To_Probabilities(100.0), pow(0.99995, 3), 1e-9);
}
