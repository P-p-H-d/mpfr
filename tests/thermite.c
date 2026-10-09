/* Test file for (physicist's) hermite polynomials.

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

#define MPFR_ORTHOGONAL_POLY_FN mpfr_hermite
#define X_LOWER_BOUND -128.0
#define X_HIGHER_BOUND 128.0
#include "torthopoly.c"

#define ARBITRARILY_LOW_PREC  10
#define RANDOM_TESTS_BATCH    20
#define DYADIC_BOUND          10
#define GENERIC_UI_RAND_MOD   10

typedef struct {
  unsigned n;
  double x;
  int sign;  /* +/-1 */
} test_case_t;

/* for n == 0:
     - H_0(Inf)  -> 1.0
     - H_0(-Inf) -> 1.0
     - H_0(NaN)  -> 1.0
   for n >= 1:
     - H_n(Inf)  -> Inf
     - H_n(-Inf) -> (-1)^n*Inf
     - H_n(NaN)  -> NaN */
static void
test_singular_input (void)
{
  mpfr_t res, x;
  int i, ret;
  unsigned n, degrees[] = { 0, 1, 2, 5, 12, 21, 30, 50 };

  mpfr_init2 (res, 200);
  mpfr_init2 (x, 200);

  for (i = 0; i < sizeof (degrees) / sizeof (unsigned); i++)
    {
      n = degrees[i];

      mpfr_set_nan (x);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (ret != 0 || !mpfr_nan_p (res))
        {
          printf ("For x = NAN, H_%u should be NAN\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      mpfr_set_inf(x, 1);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (n == 0 && (ret != 0 || mpfr_cmp_ui (res, 1) != 0))
        {
          printf ("For x = +Inf, H_0 should be +1\ngot: ");
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      if (n > 0 && (ret != 0 || !mpfr_inf_p (res) || MPFR_IS_NEG(res)))
        {
          printf ("For x = +Inf, H_%u should be +Inf\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      mpfr_set_inf(x, -1);
      ret = mpfr_hermite (res, n, x, MPFR_RNDN);
      if (n == 0 && (ret != 0 || mpfr_cmp_ui (res, 1) != 0))
        {
          printf ("For x = -Inf, H_0 should be +1\ngot: ");
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      if (n > 0 && ((n%2)==0) && (ret != 0 || !mpfr_inf_p (res) || MPFR_IS_NEG(res)))
        {
          printf ("For x = -Inf, H_%u should be +Inf\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
      /* (n%2)==1 implies n <> 0 */
      if ((n%2)==1 && (ret != 0 || !mpfr_inf_p (res) || MPFR_IS_POS(res)))
        {
          printf ("For x = -Inf, H_%u should be -Inf\ngot: ", n);
          mpfr_dump (res);
          printf ("With return value: %d\n", ret);
          exit (1);
        }
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

static void
test_zero_odd (void)
{
  int i;
  mpfr_t x, res;

  mpfr_init2 (x,IEEE754_DOUBLE_PREC);
  mpfr_init2 (res,IEEE754_DOUBLE_PREC);

  for (i = 1; i < 100; i += 2)
    {
      mpfr_set_ui (x, 0, MPFR_RNDN);
      mpfr_hermite (res, i, x, MPFR_RNDN);
      if (!MPFR_IS_ZERO (res))
        {
          printf ("H_%d(+0) should be 0; got ", i);
          mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
          printf ("\n");
          exit (1);
        }
      if ((i % 4) == 1 && MPFR_IS_NEG (res))
        {
          printf ("H_%d(+0) should be +0; got -0\n", i);
          exit (1);
        }
      if ((i % 4) == 3 && MPFR_IS_POS (res))
        {
          printf ("H_%d(+0) should be -0; got +0\n", i);
          exit (1);
        }
      mpfr_neg (x, x, MPFR_RNDN);
      mpfr_hermite (res, i, x, MPFR_RNDN);
      if (!MPFR_IS_ZERO (res))
        {
          printf ("H_%d(-0) should be 0; got ", i);
          mpfr_out_str (stdout, 10, 0, res, MPFR_RNDD);
          printf ("\n");
          exit (1);
        }
      if ((i % 4) == 1 && MPFR_IS_POS (res))
        {
          printf ("H_%d(-0) should be -0; got +0\n", i);
          exit (1);
        }
      if ((i % 4) == 3 && MPFR_IS_NEG (res))
        {
          printf ("H_%d(-0) should be +0; got -0\n", i);
          exit (1);
        }
    }

  mpfr_clear (x);
  mpfr_clear (res);
}

/* hermite(0,x) -> 1.0 */
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

  /* The first iteration should always be 1 */
  ret = mpfr_hermite (res, 0, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, one))
    {
      printf ("The first Hermite polynomial H_0 should be exactly 1.\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (one);
  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_free_cache ();
}

/* hermite(1,x) -> 2x */
static void
test_second_iteration (void)
{
  mpfr_t x, expected, res;
  int ret;

  mpfr_init2 (x, 100);
  mpfr_init2 (res, 100);
  mpfr_init2 (expected, 200);

  mpfr_set_d (x, 1.0/3.0, MPFR_RNDD);
  mpfr_mul_ui (expected, x, 2, MPFR_RNDN);

  /* The second iteration should always return 2x. Since prec(res) = prec(x),
     ret should be 0 */
  ret = mpfr_hermite (res, 1, x, MPFR_RNDN);
  if (ret != 0 || !mpfr_equal_p (res, expected))
    {
      printf ("H_1 should be 2x\ngot: ");
      mpfr_dump (res);
      printf ("With return value: %d\n", ret);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

static void
test_double_precision (void)
{
  mpfr_t x, expected, res;

  mpfr_init2 (x, IEEE754_DOUBLE_PREC);
  mpfr_init2 (res, IEEE754_DOUBLE_PREC);
  mpfr_init2 (expected, IEEE754_DOUBLE_PREC);

  /* H_3(3.49376) = 2.9924358881463502e2 with MPFR_RNDN */
  mpfr_set_d (x, 3.49376, MPFR_RNDN);
  mpfr_set_str (expected, "299.24358881463500800000000000000", 10, MPFR_RNDN);

  mpfr_hermite (res, 3, x, MPFR_RNDN);
  if (!mpfr_equal_p (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  /* H_6(-6.25) = 3.1102803906250000e6 with MPFR_RNDN */
  mpfr_set_d (x, -6.25, MPFR_RNDN);
  mpfr_set_str (expected, "3.11028039062500000000000000000e6", 10, MPFR_RNDN);

  mpfr_hermite (res, 6, x, MPFR_RNDN);
  if (!mpfr_equal_p (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  /* H_2(0.0001) = -1.9999999600000000e0 with MPFR_RNDN */
  mpfr_set_d (x, 0.0001, MPFR_RNDN);
  mpfr_set_str (expected, "-1.999999960000000000000000000000", 10, MPFR_RNDN);

  mpfr_hermite (res, 2, x, MPFR_RNDN);
  if (!mpfr_equal_p (res, expected))
    {
      printf ("test_double_precision failed;\n");
      DUMP_NUMBERS (expected, res);
      exit (1);
    }

  mpfr_clear (res);
  mpfr_clear (x);
  mpfr_clear (expected);
  mpfr_free_cache ();
}

#ifndef MPFR_USE_MINI_GMP

/* Exhaustive test of degree-n Hermite polynomial with all fractions
   a/2^b with |a| <= A and 0 <= b <= B, for precision p.
   Assume n >= 1. */
static void
test_exact (int n, int A, int B, mpfr_prec_t p)
{
  mpq_t *P0, *P1, t, u;
  int i, j, a, b, rnd;
  mpfr_t x, y, z;

  P0 = (mpq_t *) tests_allocate ((n + 1) * sizeof (mpq_t));
  P1 = (mpq_t *) tests_allocate ((n + 1) * sizeof (mpq_t));
  for (i = 0; i <= n; i++)
    {
      mpq_init (P0[i]); /* set to 0 */
      mpq_init (P1[i]); /* set to 0 */
    }
  mpq_init (t);
  mpq_init (u);
  /* use the physicist's Permite recurrence:
     H_0(x) = 1, H_1(x) = 2x,
     H_j(x) = 2x * H_{j-1}(x) - 2(j-1) * H_{j-2}(x)
     In coefficient form:
     H_j[i] = 2 * H_{j-1}[i-1] - 2*(j-1) * H_{j-2}[i] */
  mpq_set_ui (P0[0], 1, 1); /* H_0 = 1 */
  mpq_set_ui (P1[1], 2, 1); /* H_1 = 2x */
  for (j = 2; j <= n; j++)
    {
      /* P[j] = 2x * P[j-1] - 2(j-1) * P[j-2]
         thus P[j][i] = 2 * P[j-1][i-1] - 2*(j-1) * P[j-2][i].
         Invariant: P[j-2] is stored in P0, and P[j-1] in P1. */
      for (i = 0; i <= j; i++)
        {
          if (i == 0)
            mpq_set_ui (t, 0, 1);
          else
            {
              mpq_set_ui (t, 2, 1);
              mpq_mul (t, t, P1[i-1]);
            }
          /* t = 2 * P[j-1][i-1] */
          mpq_set_ui (u, 2*(j-1), 1);
          mpq_mul (u, u, P0[i]);
          /* u = 2*(j-1) * P[j-2][i] */
          mpq_sub (P0[i], t, u);
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
        RND_LOOP (rnd)\
          {
            mpfr_rnd_t r = (mpfr_rnd_t) rnd;
            mpfr_set_q (y, t, (mpfr_rnd_t) rnd); /* expected result */
            mpfr_hermite (z, n, x, r);
            if (!mpfr_equal_p (y, z))
              {
                printf ("Error in test_exact for n=%d a=%d b=%d p=%lu rnd=%s\n",
                        n, a, b, p, mpfr_print_rnd_mode (r));
                DUMP_NUMBERS (y, z);
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

#endif

/* Test overflow with reduced emax = 2.
   Test with several (n, x) pairs and all rounding modes. */
static void
test_overflow (void)
{
  mpfr_t x, y, res;
  int rnd, inex, ncases, i, is_towards_zero;
  mpfr_exp_t old_emax;
  test_case_t *curr_case;
  test_case_t cases[] =
    {
      { 2,  2.0,  1 },  /* H_2(2) = 14,     MPFR_EXP = 4  */
      { 3,  2.0,  1 },  /* H_3(2) = 40,     MPFR_EXP = 6  */
      { 2, -1.5,  1 },  /* H_2(-1.5) = 7,   MPFR_EXP = 3  */
      { 3, -1.5, -1 },  /* H_3(-1.5) = -9,  MPFR_EXP = 4  */
      { 5, -3.0, -1 },  /* H_5(-3) = -3816, MPFR_EXP = 11 */
      { 4, -2.0,  1 },  /* H_4(-2) = 76,    MPFR_EXP = 7  */
    };

  ncases = sizeof (cases) / sizeof (cases[0]);
  old_emax = mpfr_get_emax ();

  mpfr_init2 (x, 8);
  mpfr_init2 (y, 8);
  mpfr_init2 (res, 8);

  set_emax (2);

  /* y = 4 is the largest finite number: 0.11111111 * 2^emax */
  mpfr_set_ui_2exp (y, 1, 2, MPFR_RNDN);
  /* maxnum = 4 - ulp */
  mpfr_nextbelow (y);

  for (i = 0; i < ncases; i++)
    {
      curr_case = cases + i;

      RND_LOOP (rnd)
        {
          mpfr_set_d (x, curr_case->x, MPFR_RNDN);
          mpfr_clear_flags ();

          inex = mpfr_hermite (res, curr_case->n, x, (mpfr_rnd_t) rnd);

          if (!mpfr_overflow_p ())
            {
              printf ("Error in test_overflow (n = %u, x = %g, rnd = %s). "
                      "The overflow flag is not set.\n",
                      curr_case->n, curr_case->x,
                      mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
              exit (1);
            }

          /* We skip faithful rounding */
          if (rnd == MPFR_RNDF)
            continue;

          /* For rounding towards zero (RNDZ and RNDD for positive, or
             RNDZ and RNDU for negative), the result is maxnum.
             Otherwise +\-Inf (depending on the sign of x) */
          is_towards_zero = (rnd == MPFR_RNDZ);
          if (curr_case->sign > 0)
            is_towards_zero |= (rnd == MPFR_RNDD);
          else
            is_towards_zero |= (rnd == MPFR_RNDU);

          if (is_towards_zero)
            {
              if (inex * curr_case->sign >= 0)
                {
                  printf ("Error in test_overflow (n = %u, x = %g, rnd = %s). "
                          "inex has wrong sign: %d (expected_sign = %d).\n",
                          curr_case->n, curr_case->x,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd),
                          inex, -curr_case->sign);
                  exit (1);
                }
              /* Result should be maxnum with correct sign */
              if (mpfr_inf_p (res))
                {
                  printf ("Error in test_overflow (n = %u, x = %g, rnd = %s). "
                          "Got Inf, expected maxnum.\n",
                          curr_case->n, curr_case->x,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  exit (1);
                }
            }
          else
            {
              if (inex * curr_case->sign <= 0)
                {
                  printf ("Error in test_overflow (n = %u, x = %g, rnd = %s), "
                          "inex has wrong sign: %d (expected_sign = %d).\n",
                          curr_case->n, curr_case->x,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd),
                          inex, curr_case->sign);
                  exit (1);
                }
              /* Result should be +/-Inf */
              if (!mpfr_inf_p (res))
                {
                  printf ("Error in test_overflow (n = %u, x = %g, rnd = %s). Got ",
                          curr_case->n, curr_case->x,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  mpfr_dump (res);
                  printf ("Instead of %sInf.\n",
                          (curr_case->sign > 0) ? "+" : "-");
                  exit (1);
                }

              if (MPFR_SIGN (res) != curr_case->sign)
                {
                  printf ("Error in test_overflow (n = %u, x = %g, rnd = %s): "
                          "Wrong sign on Inf: ", curr_case->n, curr_case->x,
                          mpfr_print_rnd_mode ((mpfr_rnd_t) rnd));
                  exit (1);
                }
            }
        }
    }

  set_emax (old_emax);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (res);
}

/* Wrapper for tgeneric_ui.c, which calls TEST_FUNCTION(y, x, u, rnd),
   but mpfr_hermite has the integer argument before x.
   We also need to bound n because:
   - n must be non-negative (mpfr_hermite asserts n >= 0);
   - very large n would be too slow to compute for test purposes. */
static int
mpfr_hermite_generic_ui (mpfr_ptr y, mpfr_srcptr x, long n, mpfr_rnd_t rnd)
{
  n = (long) ((unsigned long) n % GENERIC_UI_RAND_MOD);
  return mpfr_hermite (y, n, x, rnd);
}

#define TEST_FUNCTION mpfr_hermite_generic_ui
#define TEST_FUNCTION_NAME "mpfr_hermite"
#define INTEGER_TYPE long
#define INT_RAND_FUNCTION() \
        (long) (randlimb () % GENERIC_UI_RAND_MOD)
#include "tgeneric_ui.c"
#undef TEST_FUNCTION
#undef TEST_FUNCTION_NAME
#undef INTEGER_TYPE

#define DEFN(N)                                                        \
  static int mpfr_hermite##N (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t r) \
  { return mpfr_hermite (y, N, x, r); }

DEFN(0)
DEFN(1)
DEFN(2)
DEFN(3)
DEFN(4)
DEFN(5)

#define TEST_FUNCTION mpfr_hermite0
#define test_generic test_generic_hermite0
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_hermite1
#define test_generic test_generic_hermite1
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_hermite2
#define test_generic test_generic_hermite2
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_hermite3
#define test_generic test_generic_hermite3
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_hermite4
#define test_generic test_generic_hermite4
#include "tgeneric.c"

#define TEST_FUNCTION mpfr_hermite5
#define test_generic test_generic_hermite5
#include "tgeneric.c"

#undef TEST_FUNCTION
#undef TEST_RANDOM_EMIN
#undef TEST_RANDOM_EMAX

int
main (void)
{
  tests_start_mpfr ();

  /* for singular input (+/-Inf or NaN), mpfr_hermite should always set res
     to NaN and return 0 */
  test_singular_input ();

  test_zero_odd ();

  /* the first two iterations are tested separately because they are the base
     cases of the recursion algorithm used to calculate mpfr_hermite */
  test_first_iteration ();
  test_second_iteration ();

  /* test physicist's Hermite polynomials with fixed test cases in double
     precision */
  test_double_precision ();

  /* perform RANDOM_TESTS_BATCH tests for each of the 10 degree chosen between
     the ones allowed  C++ standard [0, 128] */
  random_poly_suite (RANDOM_TESTS_BATCH, IEEE754_DOUBLE_PREC);

#ifndef MPFR_USE_MINI_GMP
  test_exact_dyadic ();
#endif

  /* test some specific cases of overflow in a reduced exponent range */
  test_overflow ();

  test_generic_ui (ARBITRARILY_LOW_PREC, IEEE754_DOUBLE_PREC, 6);

  test_generic_hermite0 (MPFR_PREC_MIN, 100, 10);
  test_generic_hermite1 (MPFR_PREC_MIN, 100, 10);
  test_generic_hermite2 (MPFR_PREC_MIN, 100, 100);
  test_generic_hermite3 (MPFR_PREC_MIN, 100, 100);
  test_generic_hermite4 (MPFR_PREC_MIN, 100, 100);
  test_generic_hermite5 (MPFR_PREC_MIN, 100, 100);

  tests_end_mpfr ();
  return 0;
}

#undef MPFR_ORTHOGONAL_POLY_FN
#undef X_LOWER_BOUND
#undef X_HIGHER_BOUND
#undef RANDOM_TESTS_BATCH
