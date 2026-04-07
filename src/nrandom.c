/* mpfr_nrandom (rop, state, rnd_mode) -- Generate a normal deviate with mean 0
   and variance 1 and round it to the precision of rop according to the given
   rounding mode.

Copyright 2013-2026 Free Software Foundation, Inc.
Contributed by Charles Karney <karney@alum.mit.edu>, SRI International.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.
If not, see <https://www.gnu.org/licenses/>. */

/*
 * Sampling from the normal distribution with zero mean and unit variance.
 * This uses Algorithm N given in:
 *   Charles F. F. Karney,
 *   "Sampling exactly from the normal distribution",
 *   ACM Trans. Math. Software 42(1), 3:1-14 (Jan. 2016).
 *   https://dx.doi.org/10.1145/2710016
 *   https://arxiv.org/abs/1303.6257
 *
 * mpfr_nrandom_v1 closely follows the C++ one given in the paper
 * above.  However, here, C is simplified by using gmp_urandomm_ui; the initial
 * rejection step in H just tests the leading bit of p; and the assignment of
 * the sign to the deviate using gmp_urandomb_ui.
 *
 * Improvements to this algorithm are given in:
 *   Yusong Du, Baoying Fan, and Baodian Wei,
 *   "An improved exact sampling algorithm for the standard normal
 *   distribution",
 *   Computational Statistics 37(2), 721-737 (Apr. 2022).
 *   https://doi.org/10.1007/s00180-021-01136-w
 *   https://arxiv.org/abs/2008.03855
 *
 * These improvements shave 13% off the running time for low precision results
 * (prec <= 64).  (For higher precision, the least significant bits are just
 * copied from the random number stream for both versions.)  But, as important,
 * the resulting code is simpler and easier to understand.  mpfr_nrandom_v2
 * includes these improvements as detailed in
 *   https://exrandom.sourceforge.net/3.0/algorithm.html
 *
 * There are a few "weasel words" regarding the accuracy of this
 * implementation.  The algorithm produces exactly rounded normal deviates
 * provided that gmp's random number engine delivers truly random bits.  If it
 * did, the algorithm would be perfect; however, this implementation would have
 * problems, e.g., in that the integer part of the normal deviate is
 * represented by an unsigned long, whereas in reality the integer part in
 * unbounded.  In this implementation, asserts catch overflow in the integer
 * part and similar (very, very) unlikely events.  In reality, of course, gmp's
 * random number engine has a finite internal state (19937 bits in the case of
 * the MT19937 method).  This means that these unlikely events in fact won't
 * occur.  If the asserts are triggered, then this is an indication that the
 * random number engine is defective.  (Even if a hardware random number
 * generator were used, the most likely explanation for the triggering of the
 * asserts would be that the hardware generator was broken.)
 */

#include "random_deviate.h"

/* Algorithm H: true with probability exp(-1/2).  Used by v1 and v2. */
static int
half_exp_bern (gmp_randstate_t r,
               mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  mpfr_random_deviate_reset (p);
  if (mpfr_random_deviate_tstbit (p, 1, r))
    return 1;
  for (;;)
    {
      mpfr_random_deviate_reset (q);
      if (!mpfr_random_deviate_less (q, p, r))
        return 0;
      mpfr_random_deviate_reset (p);
      if (!mpfr_random_deviate_less (p, q, r))
        return 1;
    }
}

/* Step N1: return n >= 0 with prob. exp(-n/2) * (1 - exp(-1/2)).
   Used by v1 and v2. */
static unsigned long
half_exp_geom (gmp_randstate_t r,
               mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  unsigned long n = 0;

  while (half_exp_bern (r, p, q))
    {
      ++n;
      /* Catch n wrapping around to 0; for a 32-bit unsigned long, the
       * probability of this is exp(-2^30)). */
      MPFR_ASSERTN (n != 0UL);
    }
  return n;
}

/* Step N2_v1: true with probability exp(-m*n/2).  Used by v1 only. */
static int
P (unsigned long m, unsigned long n, gmp_randstate_t r,
   mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries.  m*n is passed as two separate parameters to deal
   * with the case where m*n overflows an unsigned long.  This may be called
   * with m = 0 and n = (unsigned long)(-1) and, because m in handled in to the
   * outer loop, this routine will correctly return 1. */
  while (m--)
    {
      unsigned long k = n;
      while (k--)
        {
          if (!half_exp_bern (r, p, q))
            return 0;
        }
    }
  return 1;
}

/* Algorithm C: return (-1, 0, 1) with prob (1/m, 1/m, 1-2/m).
   Used by v1 only. */
static int
C (unsigned long m, gmp_randstate_t r)
{
  unsigned long n =  gmp_urandomm_ui (r, m);
  return n == 0 ? -1 : (n == 1 ? 0 : 1);
}

/* Algorithm B: true with prob exp(-x * (2*k + x) / (2*k + 2)).
   Used by v1 only. */
static int
B (unsigned long k, mpfr_random_deviate_t x, gmp_randstate_t r,
   mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */

  unsigned long m = 2 * k + 2;
  /* n tracks the parity of the loop; s == 1 on first trip through loop. */
  unsigned n = 0, s = 1;
  int f;

  /* Check if 2 * k + 2 would overflow; for a 32-bit unsigned long, the
   * probability of this is exp(-2^61)).  */
  MPFR_ASSERTN (k < ((unsigned long)(-1) >> 1));

  for (;; ++n, s = 0)           /* overflow of n is innocuous */
    {
      if ( ((f = k ? 0 : C (m, r)) < 0) ||
           (mpfr_random_deviate_reset (q),
            !mpfr_random_deviate_less (q, s ? x : p, r)) ||
           ((f = k ? C (m, r) : f) < 0) ||
           (f == 0 &&
            (mpfr_random_deviate_reset (p),
             !mpfr_random_deviate_less (p, x, r))) )
        break;
      mpfr_random_deviate_swap (p, q); /* an efficient way of doing p = q */
    }
  return (n & 1U) == 0;
}

/* Step N2_v2: return square root of n if it's a perfect square else -1.
   Used by v2 only.
   The probability distribution of n is exp(-n/2), so that realistically
   the maximum value of n is 100 or so, thus a naive search is enough.
   Moreover, step N1 is in Theta(n) while this step is in Theta(sqrt(n)).
   This could easily be merged into step N1 (half_exp_geom), but half_exp_geom
   is also used by v1 -- so skip this for now. */
static long
int_sqrt (unsigned long n)
{
  unsigned long k, k2, k3;

  /* k2 >= k3 test is to guard against overflow in k2 += 2*k - 1 */
  for (k = 0, k2 = 0, k3 = 0;
       k2 <= n && k2 >= k3;
       ++k, k3 = k2, k2 += 2*k - 1)
    {
      /* Here k2 = k * k; note that k^2 - (k - 1)^2 = 2*k - 1 */
      if (n == k2)
        {
          /* k = sqrt(n), thus it fits in a long. */
          MPFR_ASSERTD (k <= LONG_MAX);
          return (long) k;
        }
    }
  return -1;
}

/* Algorithm E: true with probability exp(-x) for x in (0, 1).
   This same routine appears in erandom.c, where it's called E. */
static int
trunc_exp_bern (mpfr_random_deviate_t x, gmp_randstate_t r,
                mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  mpfr_random_deviate_reset (p);
  if (!mpfr_random_deviate_less (p, x, r))
    return 1;
  for (;;)
    {
      mpfr_random_deviate_reset (q);
      if (!mpfr_random_deviate_less (q, p, r))
        return 0;
      mpfr_random_deviate_reset (p);
      if (!mpfr_random_deviate_less (p, q, r))
        return 1;
    }
}

/* Algorithm B_v2: true with prob exp(-x^2/2) for x in (0,1).  This is a
   specialization of v1's Algirithm B with k = 0.  Used by v2 only.  Here
   are the steps:

   U is a decreasing sequence x > U1 > U2 > ...
   V is a sequence V1, V2, ... with Vi < x.

   B1 [Initialize loop.] Set y = x, n = 0.
   B2 [Generate and test next samples.]
      (a) [The coin toss] With probability 1/2, go to step B4.
      (b) [The U sequence] Sample z = U; go to step B4, unless z < y.
      (c) [The V sequence] Sample r = U; go to step B4, unless r < x.
   B3 [Increment loop counter and repeat.] Set y = z, n = n+1; go to step B2.
   B4 [Test length of runs.] Set B = (n is even).

   Step B2 is entered with probability 1/2^n * x^n/n! * x^n.
*/
static int
trunc_norm_bern (mpfr_random_deviate_t x, gmp_randstate_t r,
                 mpfr_random_deviate_t p, mpfr_random_deviate_t q)
{
  /* p and q are temporaries */
  /* n tracks the parity of the loop; s == 1 on first trip through loop. */
  unsigned n = 0, s = 1;

  for (;; ++n, s = 0)           /* overflow of n is innocuous */
    {
      if (gmp_urandomb_ui (r, 1)                           /* Step B2(a) */
          || (mpfr_random_deviate_reset (q),
              !mpfr_random_deviate_less (q, s ? x : p, r)) /* Step B2(b) */
          || ((mpfr_random_deviate_reset (p),
               !mpfr_random_deviate_less (p, x, r)))       /* Step B2(c) */
          )
        break;
      mpfr_random_deviate_swap (p, q); /* an efficient way of doing p = q */
    }
  return (n & 1U) == 0;
}

/* return a normal random deviate with mean 0 and variance 1 as a MPFR.
   Version 1  */
int
mpfr_nrandom_v1 (mpfr_ptr z, gmp_randstate_t r, mpfr_rnd_t rnd)
{
  mpfr_random_deviate_t x, p, q;
  int inex;
  unsigned long k, j;

  mpfr_random_deviate_init (x);
  mpfr_random_deviate_init (p);
  mpfr_random_deviate_init (q);
  for (;;)
    {
      k = half_exp_geom (r, p, q);                   /* step 1 */
      if (!P (k, k - 1, r, p, q))
        continue;                                    /* step 2 */
      mpfr_random_deviate_reset (x);                 /* step 3 */
      for (j = 0; j <= k && B (k, x, r, p, q); ++j); /* step 4 */
      if (j > k)
        break;
    }
  mpfr_random_deviate_clear (q);
  mpfr_random_deviate_clear (p);
  /* steps 5, 6, 7 */
  inex = mpfr_random_deviate_value (gmp_urandomb_ui (r, 1), k, x, z, r, rnd);
  mpfr_random_deviate_clear (x);
  return inex;
}

/* return a normal random deviate with mean 0 and variance 1 as a MPFR.
   Version 2.  Here are the steps:

   N1 [Sample square of integer part of deviate n = k^2.]  Select integer
      n >= 0 with probability exp(−n/2) * (1 − 1/sqrt(e)).
   N2 [Testing whether n is a perfect square.] If n is a perfect square set
      k = sqrt(n); otherwise go to step N1.
   N3 [Sample fractional part of deviate x.] Set x = U.
   N4 [First adjustment of the relative probability of x.] Accept x with
      probability exp(−k*x); otherwise go to step N1.
   N5 [Second adjustment of the relative probability of x.] Accept x with
      probability exp(−x^2/2); otherwise go to step N1.
   N6 [Combine integer and fraction and assign a sign.] Set y = k + x;
      then with probability 1/2, set y = −y.
   N7 [Return result.] Set N = y.
*/
int
mpfr_nrandom_v2 (mpfr_ptr z, gmp_randstate_t r, mpfr_rnd_t rnd)
{
  mpfr_random_deviate_t x, p, q;
  int inex;
  unsigned long k;
  long j;

  mpfr_random_deviate_init (x);
  mpfr_random_deviate_init (p);
  mpfr_random_deviate_init (q);
  for (;;)
    {
      k = half_exp_geom (r, p, q);                   /* step 1 */
      if ((j = int_sqrt (k)) < 0) continue;          /* step 2 */
      k = (unsigned long) j;
      mpfr_random_deviate_reset (x);                 /* step 3 */
      while (j-- && trunc_exp_bern (x, r, p, q)) {}  /* step 4 */
      if (! (j < 0)) continue;
      if (! trunc_norm_bern (x, r, p, q) ) continue; /* step 5 */
      break;
    }
  mpfr_random_deviate_clear (q);
  mpfr_random_deviate_clear (p);
  /* steps 6, 7 */
  inex = mpfr_random_deviate_value (gmp_urandomb_ui (r, 1), k, x, z, r, rnd);
  mpfr_random_deviate_clear (x);
  return inex;
}

/* return a normal random deviate with mean 0 and variance 1 as a MPFR.
   Select the default version (currently Version 1). */
int
mpfr_nrandom (mpfr_ptr z, gmp_randstate_t r, mpfr_rnd_t rnd)
{
  return mpfr_nrandom_v1 (z, r, rnd);
}
