/* Test file for legendre polynomials.

Copyright 2025-2026 Free Software Foundation, Inc.
Contributed by Matteo Nicoli and Paul Zimmermann.

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

#include "mpfr-test.h"

#include <time.h>

#define ARBITRARILY_LOW_PREC 10
#define IEEE754_SINGLE_PREC  24
#define IEEE754_DOUBLE_PREC  53
#define MPFR_PREC_100        100
#define MPFR_PREC_200        200

#define RANDOM_TESTS_N_DEGREE 10
#define RANDOM_TESTS_BATCH    20
#define DYADIC_BOUND          35
#define GENERIC_UI_RAND_MOD   10

#define DUMP_NUMBERS(expected, got)   \
          do                          \
            {                         \
              printf ("expected: ");  \
              mpfr_dump (expected);   \
              printf ("got:      ");  \
              mpfr_dump (got);        \
            } while (0)

static const unsigned degrees[] =
{
  2,       /* The first even degree after the base cases */
  3,       /* The first odd degree after the base cases */
  10,
  50,
  128,     /* The maximum degree officially supported by C++17 */
  1024,    /* 2^10 */
};
static const char *expected_vals[] =
{
  "-0.00100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.01110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.00110000001011111100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
  "-0.11111110011011111010011011101111010001100110010000011011100111001011011101000001101100101111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-5",
  "-0.10100000000001110010100100000001010100100100100110100000000110011010001001111001010101100101000100110011011001000100100000001011011110011111100010011110111000001000111000111001110011111111101000111101e-5",
  "-0.10011011001011010000011110100001001010000100011101101001111011101111100001001000000001001001111000111010100010110101101100101110100011000100010001101010010001111010101001110011010101011000110011001001e-5",
  "-0.10100000101010101110101000101101101011000101000000110110011000110011100011010000111011000101100110010101010100011110011001001111001111110111101101010000100110001110101100110000010110100001001110000100e-8",
  "-0.10011011001100111110110010111111011011111101111001101010111101000000001111110101001101100011100101011001000110111010100000110000101101111101100001001111100111001111111011111101111001000000011000100010e-10",
};

static const int n_degrees_test = sizeof(degrees)/sizeof(unsigned);

static void
test_domain (void)
{
  int i, ret;
  mpfr_t res, upper, lower, inner, outer;

  mpfr_init2 (res, 200);

  /* init the constants to (respectively):
       - the upper bound of the domain
       - the lower bound of the domain
       - a number within the dmain
       - a number outside of the domain */
  mpfr_init2 (upper, 10);
  mpfr_set_ui (upper, 1, MPFR_RNDD);
  mpfr_init2 (lower, 200);
  mpfr_set_d (lower, -1.0, MPFR_RNDD);
  mpfr_init2 (inner, 200);
  mpfr_set_d (inner, 1.0/5.0, MPFR_RNDD);
  mpfr_init2 (outer, 200);
  /* 1.00000000001 */
  mpfr_set_d (outer, 1.0 + 1e-11, MPFR_RNDD);

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, upper, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          printf ("Upper bound input value ");
          mpfr_out_str (stdout, 10, 0, upper, MPFR_RNDD);
          printf (" should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, lower, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          printf ("Lower bound input value ");
          mpfr_out_str (stdout, 10, 0, lower, MPFR_RNDD);
          printf (" should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      mpfr_legendre (res, i, inner, MPFR_RNDN);
      if (MPFR_IS_NAN (res))
        {
          printf ("input number ");
          mpfr_out_str (stdout, 10, 0, inner, MPFR_RNDD);
          printf (" should *not* lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }
    }

  for (i = 0; i < 10; i++)
    {
      ret = mpfr_legendre (res, i, outer, MPFR_RNDN);
      if (!MPFR_IS_NAN (res))
        {
          printf ("input number outside of the domain ");
          mpfr_out_str (stdout, 10, 0, outer, MPFR_RNDD);
          printf (" should lead to a NAN result (degree: %d)\n", i);
          exit (1);
        }

      if (ret)
        {
          printf ("NaN result, should lead to an exact return value. "
                  "got: %d\n", ret);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (upper);
  mpfr_clear (lower);
  mpfr_clear (inner);
  mpfr_clear (outer);
  mpfr_free_cache ();
}

static void
test_domain_bounds (void)
{
  mpfr_t one, minus_one, res;
  int even_degree = 2, odd_degree = 3;

  mpfr_init2 (res, 200);

  /* init 1 constant */
  mpfr_init2 (one, 200);
  mpfr_set_ui (one, 1, MPFR_RNDD);

  /* init -1 constant */
  mpfr_init2 (minus_one, 200);
  mpfr_set_si (minus_one, -1, MPFR_RNDD);

  if (mpfr_legendre (res, even_degree, one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      printf ("P%d(1) should be 1; got ", even_degree);
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
      printf ("\n");
      exit (1);
    }

  if (mpfr_legendre (res, odd_degree, one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      printf ("P%d(1) should be 1; got ", odd_degree);
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
      printf ("\n");
      exit (1);
    }

  if (mpfr_legendre (res, even_degree, minus_one, MPFR_RNDD) != 0
      || !mpfr_equal_p (one, res))
    {
      printf ("P%d(-1) should be 1; got ", even_degree);

      printf ("\n");
      exit (1);
    }

  if (mpfr_legendre (res, odd_degree, minus_one, MPFR_RNDD) != 0
      || !mpfr_equal_p (minus_one, res))
    {
      printf ("P%d(-1) should be -1; got ", odd_degree);
      mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
      printf ("\n");
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (minus_one);
  mpfr_clear (res);
  mpfr_free_cache ();
}

/* legendre(0,x) -> 1.0 */
static void
test_first_iteration (void)
{
  mpfr_t one, res, x;
  int ret;

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  mpfr_set_d (x, 0.94, MPFR_RNDD);

  /* init one constant */
  mpfr_init2 (one, 200);
  mpfr_set_ui (one, 1, MPFR_RNDD);

  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  /* The first iteration should always be 1 */
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      printf ("The first legendre polynomial P_0 should be exactly 1.\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  /* The first iteration should ignore x, but still check
     that x is within the canonical domain [-1,1]*/
  mpfr_set_zero (x, 1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      printf ("P_0 should be 1, if x = +0.0\n got: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_zero (x, -1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      printf ("P_0 should be 1, if x = -0.0\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  /* For x equal to NAN, +Inf or -Inf, the result of the first iteration
     should be NAN */
  mpfr_set_nan (x);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      printf ("For x = NAN, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, 1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      printf ("For x = +Inf, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, -1);
  ret = mpfr_legendre (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_nan_p (res))
    {
      printf ("For x = -Inf, P_0 should be NAN\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* legendre(1,x) -> x */
static void
test_second_iteration (void)
{
  mpfr_t res, x;
  int ret;

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  mpfr_set_d (x, 1.0/3.0, MPFR_RNDD);

  /* The second iteration should return x itself. Since prec(res) = prec(x),
     ret should be 0 */
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, x))
    {
      printf ("P_1 should be exactly x itself\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  /* For x begin NAN, +Inf or -Inf, the result of the second iteration
     should be NAN as well */
  mpfr_set_nan (x);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (!mpfr_nan_p (res))
    {
      printf ("For x = NAN, P_1 should be NAN\n got: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  if (ret != 0)
    {
      printf ("NaN result, should lead to an exact return value. "
              "got: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, 1);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (!mpfr_nan_p (res))
    {
      printf ("For x = +Inf, P_1 should be NAN\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  if (ret != 0)
    {
      printf ("NaN result, should lead to an exact return value. "
              "got: %d\n", ret);
      exit (1);
    }
  mpfr_set_inf(x, -1);
  ret = mpfr_legendre (res, 1, x, MPFR_RNDN);
  if (!mpfr_nan_p (res))
    {
      printf ("For x = -Inf, P_1 should be NAN\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }
  if (ret != 0)
    {
      printf ("NaN result, should lead to an exact return value. "
              "got: %d\n", ret);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

static void
test_sample_with_precision (unsigned x_prec, unsigned res_prec)
{
  unsigned i;
  mpfr_t res, x, expected;
  const char *x_val = "0.5";

  mpfr_init2 (x, x_prec);
  mpfr_init2 (expected, res_prec);
  mpfr_init2 (res, res_prec);

  mpfr_set_str (x, x_val, 10, MPFR_RNDN);

  for (i = 0; i < n_degrees_test; i++)
    {
      mpfr_set_str (expected, expected_vals[i], 2, MPFR_RNDN);
      mpfr_legendre (res, degrees[i], x, MPFR_RNDN);

      if (!mpfr_eq (res, expected,  mpfr_min_prec (expected)))
        {
          printf ("Wrong value for P%d(%s) [prec(x)=%u prec(res)=%u]\n",
                  degrees[i], x_val, x_prec, res_prec);
          DUMP_NUMBERS (expected, res);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

static void
test_round (void)
{
  mpfr_t res, x;
  int ret, rnd;

  mpfr_init2 (x, IEEE754_DOUBLE_PREC);
  mpfr_init2 (res, IEEE754_SINGLE_PREC);

  mpfr_set_d(x, 1.0/3.0, MPFR_RNDD);

  /* */

  RND_LOOP_NO_RNDF (rnd)
    {
      ret = mpfr_legendre (res, 1, x, (mpfr_rnd_t) rnd);
      switch (rnd)
        {
          case MPFR_RNDN:
          case MPFR_RNDU:
          case MPFR_RNDA:
            if (ret <= 0)
              {
                printf ("For rnd=%s, P1(x), with x=1/3, mpfr_legendre should "
                        "return a positive ternary value, got: %d instead\n",
                        mpfr_print_rnd_mode ((mpfr_rnd_t) rnd), ret);
                exit (1);
              }
            break;
          case MPFR_RNDZ:
          case MPFR_RNDD:
            if (ret >= 0)
              {
                printf ("For rnd=%s, P1(x), with x=1/3, mpfr_legendre should "
                        "return a negative ternary value, got: %d instead\n",
                        mpfr_print_rnd_mode ((mpfr_rnd_t) rnd), ret);
                exit (1);
              }
            break;
          default:
            break;
        }
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* perform K random tests with degree n and precision p */
static void
test_random (int n, mpfr_prec_t p, unsigned long K)
{
  mpfr_t x, y, z, t;
  unsigned long k;
  int rnd;
  mpfr_init2 (x, p);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p + 20);
  mpfr_init2 (t, p);
  for (k = 0; k < K; k++)
    {
      mpfr_urandomb (x, RANDS); /* x is in [0,1] */
      mpfr_mul_ui (x, x, 2, MPFR_RNDN);
      mpfr_sub_ui (x, x, 1, MPFR_RNDN); /* now x is in [-1,1] */
      RND_LOOP_NO_RNDF (rnd)
        {
          mpfr_legendre (y, n, x, (mpfr_rnd_t) rnd);
          mpfr_legendre (z, n, x, MPFR_RNDN);
          if (mpfr_can_round (z, p + 20, MPFR_RNDN, (mpfr_rnd_t) rnd, p))
            {
              mpfr_set (t, z, (mpfr_rnd_t) rnd);
              if (mpfr_cmp (y, t))
                {
                  printf ("Error in mpfr_legendre for n=%d x=", n);
                  mpfr_out_str (stdout, 16, 0, x, MPFR_RNDN);
                  printf (" rnd=%s\n", mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  DUMP_NUMBERS (t, y);
                  exit (1);
                }
            }
        }
    }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

static void
random_array (int *array, int size, int inf, int sup)
{
  int i, tmp;
  for (i = 0; i < size; i++)
    {
      if (inf > sup)
        {
          tmp = inf;
          inf = sup;
          sup = tmp;
        }
      array[i] = inf + rand() % (sup - inf + 1);
    }
}

static void
random_test_suite (int num_degrees, int num_tests)
{
  /* we set the minimum degree to 2 to skip the two base cases P0 and P1,
     and the maximum degree to 128 to limit the range of degrees tested
     to the same limit of the C++ standard */
  int i, min_degree = 2, max_degree = 128;
  int *test_degrees;

  srand(time(NULL));

  test_degrees = (int *) tests_allocate (num_degrees * sizeof(int));
  if (!test_degrees)
    {
      printf ("Could not allocate memory for random tests\n");
      exit (1);
    }

  random_array (test_degrees, num_degrees, min_degree, max_degree);

  for (i = 0; i < num_degrees; i++)
    {
      test_random (test_degrees[i], IEEE754_DOUBLE_PREC, num_tests);
    }
  tests_free (test_degrees, num_degrees * sizeof(int));
}

static void
bug20251001 (void)
{
  mpfr_t x, y;
  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);

  /* bug found with GMP_CHECK_RANDOMIZE=1759307919467569:
     Error in mpfr_legendre for n=2 x=-4.7184833984611241e-1 rnd=MPFR_RNDD
     expected: -0.10101010000001100000110110100001000111001110111100001E-2
     got:      -0.10101010000001100000110110100001000111001110111100000E-2 */
  mpfr_set_str (x, "-4.7184833984611241e-1", 10, MPFR_RNDN);
  mpfr_legendre (y, 2, x, MPFR_RNDD);
  mpfr_set_str_binary (x, "-0.10101010000001100000110110100001000111001110111100001E-2");
  if (mpfr_cmp (y, x) != 0) {
    printf ("Error in bug20251001 (1)\n");
    DUMP_NUMBERS (x,y);
    exit (1);
  }

  /* bug found with GMP_CHECK_RANDOMIZE=1759308817539561:
     Error in mpfr_legendre for n=10 x=4.608a667782dd0@-1 rnd=MPFR_RNDU
     expected: 0.11111010001111111110111110101011110010000001110000001E-2
     got:      0.11111010001111111110111110101011110010000001110000000E-2 */
  mpfr_set_str (x, "4.608a667782dd0@-1", 16, MPFR_RNDN);
  mpfr_legendre (y, 10, x, MPFR_RNDU);
  mpfr_set_str_binary (x, "0.11111010001111111110111110101011110010000001110000001E-2");
  if (mpfr_cmp (y, x) != 0) {
    printf ("Error in bug20251001 (2)\n");
    DUMP_NUMBERS (x,y);
    exit (1);
  }

  /* bug found with GMP_CHECK_RANDOMIZE=1759309465303933:
     Error in mpfr_legendre for n=3 x=5.6a59d5d01d9d0@-1 rnd=MPFR_RNDZ
     expected: -0.11010010010011110110100000101001000001111000100101111E-1
     got:      -0.11010010010011110110100000101001000001111000100110000E-1 */
  mpfr_set_str (x, "5.6a59d5d01d9d0@-1", 16, MPFR_RNDN);
  mpfr_legendre (y, 3, x, MPFR_RNDZ);
  mpfr_set_str_binary (x, "-0.11010010010011110110100000101001000001111000100101111E-1");
  if (mpfr_cmp (y, x) != 0) {
    printf ("Error in bug20251001 (3)\n");
    DUMP_NUMBERS (x,y);
    exit (1);
  }

  mpfr_clear (x);
  mpfr_clear (y);
}

static void
bug20260205 (void)
{
  mpfr_t x, y;
  mpfr_init2 (x, 53);
  mpfr_init2 (y, 53);

  /* Error in mpfr_legendre for n=2 x=-9.3dbaa2fd514c0@-1 rnd=MPFR_RNDZ
     expected: 0.11001000000101110011101010100100001111111000011100111E-11
     got:      0.11001000000101110011101010100100001111111000011101000E-11 */
  mpfr_set_str (x, "-9.3dbaa2fd514c0@-1", 16, MPFR_RNDN);
  mpfr_legendre (y, 2, x, MPFR_RNDZ);
  mpfr_set_str_binary (x, "0.11001000000101110011101010100100001111111000011100111E-11");
  if (mpfr_cmp (y, x) != 0) {
    printf ("Error in bug20260205\n");
    DUMP_NUMBERS (x,y);
    exit (1);
  }

  mpfr_clear (x);
  mpfr_clear (y);
}

/* Exhaustive test of degree-n Legendre polynomial with all fractions
   a/2^b with |a| <= A and |b| <= B, for precision p.
   Assume n >= 1. */
static void
test_exact (int n, int A, int B, mpfr_prec_t p)
{
  mpq_t *P0, *P1, t, u;
  int i, j, a, b, rnd;
  mpfr_t x, y, z;

  P0 = (mpq_t*) tests_allocate ((n + 1) * sizeof (mpq_t));
  P1 = (mpq_t*) tests_allocate ((n + 1) * sizeof (mpq_t));
  for (i = 0; i <= n; i++)
    {
      mpq_init (P0[i]); /* set to 0 */
      mpq_init (P1[i]); /* set to 0 */
    }
  mpq_init (t);
  mpq_init (u);
  /* use Bonnet's recursion */
  mpq_set_ui (P0[0], 1, 1); /* P[0] = 1 */
  mpq_set_ui (P1[1], 1, 1); /* P[1] = x */
  for (j = 2; j <= n; j++)
    {
      /* P[j] = ((2j-1)*x*P[j-1]-(j-1)*P[j-2])/j
         thus P[j][i] = ((2j-1)*P[j-1][i-1] - (j-1)*P[j-2][i])/j.
         Invariant: P[j-2] is stored in P0, and P[j-1] in P1. */
      for (i = 0; i <= j; i++)
        {
          if (i == 0)
            mpq_set_ui (t, 0, 1);
          else
            {
              mpq_set_ui (t, 2*j-1, 1);
              mpq_mul (t, t, P1[i-1]);
            }
          /* t = (2j-1)*P[j-1][i-1] */
          mpq_set_ui (u, j-1, 1);
          mpq_mul (u, u, P0[i]);
          /* u = (j-1)*P[j-2][i] */
          mpq_sub (t, t, u);
          /* t = (2j-1)*P[j-1][i-1] - (j-1)*P[j-2][i] */
          mpq_set_ui (u, j, 1);
          mpq_div (P0[i], t, u);
          /* now P0[i] contains P[j][i] */
        }
      /* swap P0 and P1 */
      for (i = 0; i <= j; i++)
        mpq_swap (P0[i], P1[i]);
    }

  mpfr_init2 (x, 64);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p);

  for (a = -A; a <= A; a++)
    for (b = 0; b <= B; b++)
      {
        /* compute t = Pn(a/2^b) */
        mpq_set_si (u, a, 1ul<<b);
        mpq_set (t, P1[n]);
        for (i = n-1; i >= 0; i--)
          {
            mpq_mul (t, t, u);
            mpq_add (t, t, P1[i]);
          }

        /* now t = Pn(a/2^b) exactly */

        mpfr_set_si_2exp (x, a, -b, MPFR_RNDN);
        RND_LOOP (rnd)
          {
            mpfr_rnd_t r = (mpfr_rnd_t) rnd;
            mpfr_set_q (y, t, (mpfr_rnd_t) rnd); /* expected result */
            mpfr_legendre (z, n, x, r);
            if (mpfr_cmp (y, z))
              {
                printf ("Error in test_exact for n=%d a=%d b=%d p=%lu rnd=%s\n",
                        n, a, b, (unsigned long) p, mpfr_print_rnd_mode (r));
                printf ("expected ");
                mpfr_out_str (stdout, 16, 0, y, MPFR_RNDN);
                printf ("\ngot      ");
                mpfr_out_str (stdout, 16, 0, z, MPFR_RNDN);
                printf ("\n");
                exit (1);
            }
          }
      }

  for (i = 0; i <= n; i++)
    {
      mpq_clear (P0[i]);
      mpq_clear (P1[i]);
    }
  tests_free (P0, (n + 1) * sizeof (mpq_t));
  tests_free (P1, (n + 1) * sizeof (mpq_t));
  mpq_clear (t);
  mpq_clear (u);
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
test_exact_dyadic (void)
{
  int n;
  mpfr_prec_t p;

  for (n = 1; n <= 10; n++)
    for (p = DYADIC_BOUND - 3; p <= DYADIC_BOUND; p++)
      test_exact (n, DYADIC_BOUND, DYADIC_BOUND, p);
}

/* test x=0 for n even */
static void
test_zero_even (void)
{
  long n;
  mpfr_prec_t p;
  mpq_t q, r;
  mpq_init (q);
  mpq_init (r);
  for (n = 0; n <= 100; n+=2)
    {
      /* invariant: q = Pn(0) for n even */
      if (n == 0)
        mpq_set_ui (q, 1, 1); /* P_0(0) = 1 */
      else
        {
          mpq_set_si (r, -(n-1), n);
          mpq_mul (q, q, r);
        }
      for (p = 1; p <= 100; p++)
        {
          mpfr_t x, y, z;
          int ret, ret2, rnd;
          mpfr_init2 (x, p);
          mpfr_init2 (y, p);
          mpfr_init2 (z, p);
          mpfr_set_ui (x, 0, MPFR_RNDN);
          RND_LOOP_NO_RNDF (rnd)
            {
              ret = mpfr_legendre (y, n, x, (mpfr_rnd_t) rnd);
              ret2 = mpfr_set_q (z, q, (mpfr_rnd_t) rnd);
              if (mpfr_cmp (y, z))
                {
                  printf ("test_zero failed for n=%ld p=%lu rnd=%s\n",
                          n, (unsigned long) p,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  DUMP_NUMBERS(z, y);
                  exit (1);
                }
              if (ret * ret2 < 0 || (!!ret ^ !!ret2))
                {
                  printf ("test_zero: incompatible ternary values for "
                          "n=%ld p=%lu rnd=%s\n",
                          n, (unsigned long) p,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  printf ("got %d expected %d\n", ret, ret2);
                  exit (1);
                }
            }
          mpfr_clear (x);
          mpfr_clear (y);
          mpfr_clear (z);
        }
    }
  mpq_clear (q);
  mpq_clear (r);
}

static void
test_zero_odd (void)
{
  int i;
  mpfr_t x, res;

  mpfr_init2 (x,IEEE754_DOUBLE_PREC);
  mpfr_init2 (res,IEEE754_DOUBLE_PREC);
  mpfr_set_ui (x, 0, MPFR_RNDN);

  for (i = 1; i < 100; i += 2)
    {
      mpfr_legendre (res, i, x, MPFR_RNDN);
      if (!MPFR_IS_ZERO (res))
        {
          printf ("P_%d(0) should be 0; got ", i);
          mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
          printf ("\n");
          exit (1);
        }

      if ((i % 4) == 1 && MPFR_IS_NEG (res))
        {
          printf ("P_%d(0) should be +0; got -0\n", i);
          exit (1);
        }

      if ((i % 4) == 3 && MPFR_IS_POS (res))
        {
          printf ("P_%d(0) should be -0; got +0\n", i);
          exit (1);
        }
    }

  mpfr_clear (x);
  mpfr_clear (res);
}

/* Wrapper for tgeneric_ui.c, which calls TEST_FUNCTION(y, x, u, rnd),
   but mpfr_legendre has the integer argument before x.
   We also need to bound n because:
   - n must be non-negative (mpfr_legendre asserts n >= 0);
   - very large n would be too slow to compute for test purposes. */
static int
mpfr_legendre_generic_ui (mpfr_ptr y, mpfr_srcptr x, long n, mpfr_rnd_t rnd)
{
  n = (long) ((unsigned long) n % GENERIC_UI_RAND_MOD);
  return mpfr_legendre (y, n, x, rnd);
}

#define TEST_FUNCTION mpfr_legendre_generic_ui
#define TEST_FUNCTION_NAME "mpfr_legendre"
#define INTEGER_TYPE long
#define INT_RAND_FUNCTION() \
        (long) (randlimb () % GENERIC_UI_RAND_MOD)
#include "tgeneric_ui.c"

int
main (void)
{
  tests_start_mpfr ();

  test_zero_even ();
  test_zero_odd ();

  /* The canonical domain of Legendre polynomials is [-1,1]. mpfr_legendre
     should return false for any x outside of the canonical domain */
  test_domain ();

  /* Upper and lower bounds are computed without using Bonnet's recursion */
  test_domain_bounds ();

  /* The first two iterations are tested separately because they are the base
     cases of the Bonnet's recursion algorithm used by mpfr_legendre */
  test_first_iteration ();
  test_second_iteration ();

  /* This test checks the ternary value returned by mpfr_legendre */
  test_round ();

  /* res precision is arbitrarily low. ARBITRARILY_LOW_PREC should be lower
     than any other precision */
  test_sample_with_precision (IEEE754_SINGLE_PREC, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (MPFR_PREC_100, ARBITRARILY_LOW_PREC);
  test_sample_with_precision (MPFR_PREC_200, ARBITRARILY_LOW_PREC);

  /* res precision is IEEE754_DOUBLE_PREC. It's higher than
     IEEE754_SINGLE_PREC and lower than MPFR_PREC_100 and
     MPFR_PREC_200 */
  test_sample_with_precision (IEEE754_SINGLE_PREC, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (MPFR_PREC_100, IEEE754_DOUBLE_PREC);
  test_sample_with_precision (MPFR_PREC_200, IEEE754_DOUBLE_PREC);

  /* res precision is MPFR_PREC_200, the highest available.
     NOTE: to reach higher precisions, the expected values should be
     recalculated, because they all contains a 200-bit significand */
  test_sample_with_precision (IEEE754_SINGLE_PREC, MPFR_PREC_200);
  test_sample_with_precision (IEEE754_DOUBLE_PREC, MPFR_PREC_200);
  test_sample_with_precision (MPFR_PREC_100, MPFR_PREC_200);
  test_sample_with_precision (MPFR_PREC_200, MPFR_PREC_200);

  /* perform RANDOM_TESTS_BATCH tests for each of the RANDOM_TESTS_N_DEGREE */
  random_test_suite (RANDOM_TESTS_N_DEGREE, RANDOM_TESTS_BATCH);

  bug20251001 ();

  bug20260205 ();

  test_exact_dyadic ();

  test_generic_ui (ARBITRARILY_LOW_PREC, IEEE754_DOUBLE_PREC, 6);

  tests_end_mpfr ();
  return 0;
}
