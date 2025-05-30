/* Test file for mpfr_sqrt.

Copyright 1999, 2001-2025 Free Software Foundation, Inc.
Contributed by the Pascaline and Caramba projects, INRIA.

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

#ifdef CHECK_EXTERNAL
static int
test_sqrt (mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd_mode)
{
  int res;
  int ok = rnd_mode == MPFR_RNDN && mpfr_number_p (b);
  if (ok)
    {
      mpfr_print_raw (b);
    }
  res = mpfr_sqrt (a, b, rnd_mode);
  if (ok)
    {
      printf (" ");
      mpfr_print_raw (a);
      printf ("\n");
    }
  return res;
}
#else
#define test_sqrt mpfr_sqrt
#endif

static void
check3 (const char *as, mpfr_rnd_t rnd_mode, const char *qs)
{
  mpfr_t q;

  mpfr_init2 (q, 53);
  mpfr_set_str1 (q, as);
  test_sqrt (q, q, rnd_mode);
  if (mpfr_cmp_str1 (q, qs) )
    {
      printf ("mpfr_sqrt failed for a=%s, rnd_mode=%s\n",
              as, mpfr_print_rnd_mode (rnd_mode));
      printf ("expected sqrt is %s, got ",qs);
      mpfr_out_str (stdout, 10, 0, q, MPFR_RNDN);
      putchar ('\n');
      exit (1);
    }
  mpfr_clear (q);
}

static void
check4 (const char *as, mpfr_rnd_t rnd_mode, const char *qs)
{
  mpfr_t q;

  mpfr_init2 (q, 53);
  mpfr_set_str1 (q, as);
  test_sqrt (q, q, rnd_mode);
  if (mpfr_cmp_str (q, qs, 16, MPFR_RNDN))
    {
      printf ("mpfr_sqrt failed for a=%s, rnd_mode=%s\n",
              as, mpfr_print_rnd_mode (rnd_mode));
      printf ("expected %s\ngot      ", qs);
      mpfr_out_str (stdout, 16, 0, q, MPFR_RNDN);
      printf ("\n");
      mpfr_clear (q);
      exit (1);
    }
  mpfr_clear (q);
}

static void
check24 (const char *as, mpfr_rnd_t rnd_mode, const char *qs)
{
  mpfr_t q;

  mpfr_init2 (q, 24);
  mpfr_set_str1 (q, as);
  test_sqrt (q, q, rnd_mode);
  if (mpfr_cmp_str1 (q, qs))
    {
      printf ("mpfr_sqrt failed for a=%s, prec=24, rnd_mode=%s\n",
              as, mpfr_print_rnd_mode(rnd_mode));
      printf ("expected sqrt is %s, got ",qs);
      mpfr_out_str (stdout, 10, 0, q, MPFR_RNDN);
      printf ("\n");
      exit (1);
    }
  mpfr_clear (q);
}

static void
check_diverse (const char *as, mpfr_prec_t p, const char *qs)
{
  mpfr_t q;

  mpfr_init2 (q, p);
  mpfr_set_str1 (q, as);
  test_sqrt (q, q, MPFR_RNDN);
  if (mpfr_cmp_str1 (q, qs))
    {
      printf ("mpfr_sqrt failed for a=%s, prec=%lu, rnd_mode=%s\n",
              as, (unsigned long) p, mpfr_print_rnd_mode (MPFR_RNDN));
      printf ("expected sqrt is %s, got ", qs);
      mpfr_out_str (stdout, 10, 0, q, MPFR_RNDN);
      printf ("\n");
      exit (1);
    }
  mpfr_clear (q);
}

/* the following examples come from the paper "Number-theoretic Test
   Generation for Directed Rounding" from Michael Parks, Table 3 */
static void
check_float (void)
{
  check24("70368760954880.0", MPFR_RNDN, "8.388609e6");
  check24("281474943156224.0", MPFR_RNDN, "1.6777215e7");
  check24("70368777732096.0", MPFR_RNDN, "8.388610e6");
  check24("281474909601792.0", MPFR_RNDN, "1.6777214e7");
  check24("100216216748032.0", MPFR_RNDN, "1.0010805e7");
  check24("120137273311232.0", MPFR_RNDN, "1.0960715e7");
  check24("229674600890368.0", MPFR_RNDN, "1.5155019e7");
  check24("70368794509312.0", MPFR_RNDN, "8.388611e6");
  check24("281474876047360.0", MPFR_RNDN, "1.6777213e7");
  check24("91214552498176.0", MPFR_RNDN, "9.550631e6");

  check24("70368760954880.0", MPFR_RNDZ, "8.388608e6");
  check24("281474943156224.0", MPFR_RNDZ, "1.6777214e7");
  check24("70368777732096.0", MPFR_RNDZ, "8.388609e6");
  check24("281474909601792.0", MPFR_RNDZ, "1.6777213e7");
  check24("100216216748032.0", MPFR_RNDZ, "1.0010805e7");
  check24("120137273311232.0", MPFR_RNDZ, "1.0960715e7");
  check24("229674600890368.0", MPFR_RNDZ, "1.5155019e7");
  check24("70368794509312.0", MPFR_RNDZ, "8.38861e6");
  check24("281474876047360.0", MPFR_RNDZ, "1.6777212e7");
  check24("91214552498176.0", MPFR_RNDZ, "9.550631e6");

  check24("70368760954880.0", MPFR_RNDU, "8.388609e6");
  check24("281474943156224.0",MPFR_RNDU, "1.6777215e7");
  check24("70368777732096.0", MPFR_RNDU, "8.388610e6");
  check24("281474909601792.0", MPFR_RNDU, "1.6777214e7");
  check24("100216216748032.0", MPFR_RNDU, "1.0010806e7");
  check24("120137273311232.0", MPFR_RNDU, "1.0960716e7");
  check24("229674600890368.0", MPFR_RNDU, "1.515502e7");
  check24("70368794509312.0", MPFR_RNDU, "8.388611e6");
  check24("281474876047360.0", MPFR_RNDU, "1.6777213e7");
  check24("91214552498176.0", MPFR_RNDU, "9.550632e6");

  check24("70368760954880.0", MPFR_RNDD, "8.388608e6");
  check24("281474943156224.0", MPFR_RNDD, "1.6777214e7");
  check24("70368777732096.0", MPFR_RNDD, "8.388609e6");
  check24("281474909601792.0", MPFR_RNDD, "1.6777213e7");
  check24("100216216748032.0", MPFR_RNDD, "1.0010805e7");
  check24("120137273311232.0", MPFR_RNDD, "1.0960715e7");
  check24("229674600890368.0", MPFR_RNDD, "1.5155019e7");
  check24("70368794509312.0", MPFR_RNDD, "8.38861e6");
  check24("281474876047360.0", MPFR_RNDD, "1.6777212e7");
  check24("91214552498176.0", MPFR_RNDD, "9.550631e6");

  /* check that rounding away is just rounding toward positive infinity */
  check24("91214552498176.0", MPFR_RNDA, "9.550632e6");
}

static void
special (void)
{
  mpfr_t x, y, z;
  int inexact;
  mpfr_prec_t p;

  mpfr_init (x);
  mpfr_init (y);
  mpfr_init (z);

  mpfr_set_prec (x, 64);
  mpfr_set_str_binary (x, "1010000010100011011001010101010010001100001101011101110001011001E-1");
  mpfr_set_prec (y, 32);
  test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 2405743844UL))
    {
      printf ("Error for n^2+n+1/2 with n=2405743843\n");
      exit (1);
    }

  mpfr_set_prec (x, 65);
  mpfr_set_str_binary (x, "10100000101000110110010101010100100011000011010111011100010110001E-2");
  mpfr_set_prec (y, 32);
  test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 2405743844UL))
    {
      printf ("Error for n^2+n+1/4 with n=2405743843\n");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_set_prec (x, 66);
  mpfr_set_str_binary (x, "101000001010001101100101010101001000110000110101110111000101100011E-3");
  mpfr_set_prec (y, 32);
  test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 2405743844UL))
    {
      printf ("Error for n^2+n+1/4+1/8 with n=2405743843\n");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_set_prec (x, 66);
  mpfr_set_str_binary (x, "101000001010001101100101010101001000110000110101110111000101100001E-3");
  mpfr_set_prec (y, 32);
  test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 2405743843UL))
    {
      printf ("Error for n^2+n+1/8 with n=2405743843\n");
      mpfr_dump (y);
      exit (1);
    }

  mpfr_set_prec (x, 27);
  mpfr_set_str_binary (x, "0.110100111010101000010001011");
  if ((inexact = test_sqrt (x, x, MPFR_RNDZ)) >= 0)
    {
      printf ("Wrong inexact flag: expected -1, got %d\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 2);
  for (p=2; p<1000; p++)
    {
      mpfr_set_prec (z, p);
      mpfr_set_ui (z, 1, MPFR_RNDN);
      mpfr_nexttoinf (z);
      test_sqrt (x, z, MPFR_RNDU);
      if (mpfr_cmp_ui_2exp(x, 3, -1))
        {
          printf ("Error: sqrt(1+ulp(1), up) should give 1.5 (prec=%u)\n",
                  (unsigned int) p);
          printf ("got "); mpfr_dump (x);
          exit (1);
        }
    }

  /* check inexact flag */
  mpfr_set_prec (x, 5);
  mpfr_set_str_binary (x, "1.1001E-2");
  if ((inexact = test_sqrt (x, x, MPFR_RNDN)))
    {
      printf ("Wrong inexact flag: expected 0, got %d\n", inexact);
      exit (1);
    }

  mpfr_set_prec (x, 2);
  mpfr_set_prec (z, 2);

  /* checks the sign is correctly set */
  mpfr_set_si (x, 1,  MPFR_RNDN);
  mpfr_set_si (z, -1, MPFR_RNDN);
  test_sqrt (z, x, MPFR_RNDN);
  if (mpfr_cmp_ui (z, 0) < 0)
    {
      printf ("Error: square root of 1 gives ");
      mpfr_dump (z);
      exit (1);
    }

  mpfr_set_prec (x, 192);
  mpfr_set_prec (z, 160);
  mpfr_set_str_binary (z, "0.1011010100000100100100100110011001011100100100000011000111011001011101101101110000110100001000100001100001011000E1");
  mpfr_set_prec (x, 160);
  test_sqrt(x, z, MPFR_RNDN);
  test_sqrt(z, x, MPFR_RNDN);

  mpfr_set_prec (x, 53);
  mpfr_set_str (x, "8093416094703476.0", 10, MPFR_RNDN);
  mpfr_div_2ui (x, x, 1075, MPFR_RNDN);
  test_sqrt (x, x, MPFR_RNDN);
  mpfr_set_str (z, "1e55596835b5ef@-141", 16, MPFR_RNDN);
  if (mpfr_cmp (x, z))
    {
      printf ("Error: square root of 8093416094703476*2^(-1075)\n");
      printf ("expected "); mpfr_dump (z);
      printf ("got      "); mpfr_dump (x);
      exit (1);
    }

  mpfr_set_prec (x, 33);
  mpfr_set_str_binary (x, "0.111011011011110001100111111001000e-10");
  mpfr_set_prec (z, 157);
  inexact = test_sqrt (z, x, MPFR_RNDN);
  mpfr_set_prec (x, 157);
  mpfr_set_str_binary (x, "0.11110110101100101111001011100011100011100001101010111011010000100111011000111110100001001011110011111100101110010110010110011001011011010110010000011001101E-5");
  if (mpfr_cmp (x, z))
    {
      printf ("Error: square root (1)\n");
      exit (1);
    }
  if (inexact <= 0)
    {
      printf ("Error: wrong inexact flag (1)\n");
      exit (1);
    }

  /* case prec(result) << prec(input) */
  mpfr_set_prec (z, 2);
  for (p = mpfr_get_prec (z); p < 1000; p++)
    {
      mpfr_set_prec (x, p);
      mpfr_set_ui (x, 1, MPFR_RNDN);
      mpfr_nextabove (x);
      /* 1.0 < x <= 1.5 thus 1 < sqrt(x) <= 1.23 */
      inexact = test_sqrt (z, x, MPFR_RNDN);
      MPFR_ASSERTN(inexact < 0 && mpfr_cmp_ui (z, 1) == 0);
      inexact = test_sqrt (z, x, MPFR_RNDZ);
      MPFR_ASSERTN(inexact < 0 && mpfr_cmp_ui (z, 1) == 0);
      inexact = test_sqrt (z, x, MPFR_RNDU);
      MPFR_ASSERTN(inexact > 0 && mpfr_cmp_ui_2exp (z, 3, -1) == 0);
      inexact = test_sqrt (z, x, MPFR_RNDD);
      MPFR_ASSERTN(inexact < 0 && mpfr_cmp_ui (z, 1) == 0);
      inexact = test_sqrt (z, x, MPFR_RNDA);
      MPFR_ASSERTN(inexact > 0 && mpfr_cmp_ui_2exp (z, 3, -1) == 0);
    }

  /* corner case rw = 0 in rounding to nearest */
  mpfr_set_prec (z, GMP_NUMB_BITS - 1);
  mpfr_set_prec (y, GMP_NUMB_BITS - 1);
  mpfr_set_ui (y, 1, MPFR_RNDN);
  mpfr_mul_2ui (y, y, GMP_NUMB_BITS - 1, MPFR_RNDN);
  mpfr_nextabove (y);
  for (p = 2 * GMP_NUMB_BITS - 1; p <= 1000; p++)
    {
      mpfr_set_prec (x, p);
      mpfr_set_ui (x, 1, MPFR_RNDN);
      mpfr_set_exp (x, GMP_NUMB_BITS);
      mpfr_add_ui (x, x, 1, MPFR_RNDN);
      /* now x = 2^(GMP_NUMB_BITS - 1) + 1 (GMP_NUMB_BITS bits) */
      inexact = mpfr_sqr (x, x, MPFR_RNDN);
      MPFR_ASSERTN (inexact == 0); /* exact */
      inexact = test_sqrt (z, x, MPFR_RNDN);
      /* even rule: z should be 2^(GMP_NUMB_BITS - 1) */
      MPFR_ASSERTN (inexact < 0);
      inexact = mpfr_cmp_ui_2exp (z, 1, GMP_NUMB_BITS - 1);
      MPFR_ASSERTN (inexact == 0);
      mpfr_nextbelow (x);
      /* now x is just below [2^(GMP_NUMB_BITS - 1) + 1]^2 */
      inexact = test_sqrt (z, x, MPFR_RNDN);
      MPFR_ASSERTN(inexact < 0 &&
                   mpfr_cmp_ui_2exp (z, 1, GMP_NUMB_BITS - 1) == 0);
      mpfr_nextabove (x);
      mpfr_nextabove (x);
      /* now x is just above [2^(GMP_NUMB_BITS - 1) + 1]^2 */
      inexact = test_sqrt (z, x, MPFR_RNDN);
      if (mpfr_cmp (z, y))
        {
          printf ("Error for sqrt(x) in rounding to nearest\n");
          printf ("x="); mpfr_dump (x);
          printf ("Expected "); mpfr_dump (y);
          printf ("Got      "); mpfr_dump (z);
          exit (1);
        }
      if (inexact <= 0)
        {
          printf ("Wrong inexact flag in corner case for p = %lu\n", (unsigned long) p);
          exit (1);
        }
    }

  mpfr_set_prec (x, 1000);
  mpfr_set_ui (x, 9, MPFR_RNDN);
  mpfr_set_prec (y, 10);
  inexact = test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 3) || inexact != 0)
    {
      printf ("Error in sqrt(9:1000) for prec=10\n");
      exit (1);
    }
  mpfr_set_prec (y, GMP_NUMB_BITS);
  mpfr_nextabove (x);
  inexact = test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 3) || inexact >= 0)
    {
      printf ("Error in sqrt(9:1000) for prec=%d\n", (int) GMP_NUMB_BITS);
      exit (1);
    }
  mpfr_set_prec (x, 2 * GMP_NUMB_BITS);
  mpfr_set_prec (y, GMP_NUMB_BITS);
  mpfr_set_ui (y, 1, MPFR_RNDN);
  mpfr_nextabove (y);
  mpfr_set (x, y, MPFR_RNDN);
  inexact = test_sqrt (y, x, MPFR_RNDN);
  if (mpfr_cmp_ui (y, 1) || inexact >= 0)
    {
      printf ("Error in sqrt(1) for prec=%d\n", (int) GMP_NUMB_BITS);
      mpfr_dump (y);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
check_inexact (mpfr_prec_t p)
{
  mpfr_t x, y, z;
  mpfr_rnd_t rnd;
  int inexact, sign;

  mpfr_init2 (x, p);
  mpfr_init2 (y, p);
  mpfr_init2 (z, 2*p);
  mpfr_urandomb (x, RANDS);
  rnd = RND_RAND_NO_RNDF ();
  inexact = test_sqrt (y, x, rnd);
  if (mpfr_sqr (z, y, rnd)) /* exact since prec(z) = 2*prec(y) */
    {
      printf ("Error: multiplication should be exact\n");
      exit (1);
    }
  mpfr_sub (z, z, x, rnd); /* exact also */
  sign = mpfr_cmp_ui (z, 0);
  if (((inexact == 0) && (sign)) ||
      ((inexact > 0) && (sign <= 0)) ||
      ((inexact < 0) && (sign >= 0)))
    {
      printf ("Error with rnd=%s: wrong ternary value, expected %d, got %d\n",
              mpfr_print_rnd_mode (rnd), sign, inexact);
      printf ("x=");
      mpfr_dump (x);
      printf ("y=");
      mpfr_dump (y);
      exit (1);
    }
  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
}

static void
check_singular (void)
{
  mpfr_t  x, got;

  mpfr_init2 (x, 100L);
  mpfr_init2 (got, 100L);

  /* sqrt(NaN) == NaN */
  MPFR_SET_NAN (x);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_nan_p (got));

  /* sqrt(-1) == NaN */
  mpfr_set_si (x, -1L, MPFR_RNDZ);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_nan_p (got));

  /* sqrt(+inf) == +inf */
  MPFR_SET_INF (x);
  MPFR_SET_POS (x);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_inf_p (got));

  /* sqrt(-inf) == NaN */
  MPFR_SET_INF (x);
  MPFR_SET_NEG (x);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_nan_p (got));

  /* sqrt(+0) == +0 */
  mpfr_set_si (x, 0L, MPFR_RNDZ);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_number_p (got));
  MPFR_ASSERTN (mpfr_cmp_ui (got, 0L) == 0);
  MPFR_ASSERTN (MPFR_IS_POS (got));

  /* sqrt(-0) == -0 */
  mpfr_set_si (x, 0L, MPFR_RNDZ);
  MPFR_SET_NEG (x);
  MPFR_ASSERTN (test_sqrt (got, x, MPFR_RNDZ) == 0); /* exact */
  MPFR_ASSERTN (mpfr_number_p (got));
  MPFR_ASSERTN (mpfr_cmp_ui (got, 0L) == 0);
  MPFR_ASSERTN (MPFR_IS_NEG (got));

  mpfr_clear (x);
  mpfr_clear (got);
}

/* check that -1 <= x/sqrt(x^2+s*y^2) <= 1 for rounding to nearest or up
   with s = 0 and s = 1 */
static void
test_property1 (mpfr_prec_t p, mpfr_rnd_t r, int s)
{
  mpfr_t x, y, z, t;

  mpfr_init2 (x, p);
  mpfr_init2 (y, p);
  mpfr_init2 (z, p);
  mpfr_init2 (t, p);

  mpfr_urandomb (x, RANDS);
  mpfr_mul (z, x, x, r);
  if (s)
    {
      mpfr_urandomb (y, RANDS);
      mpfr_mul (t, y, y, r);
      mpfr_add (z, z, t, r);
    }
  mpfr_sqrt (z, z, r);
  mpfr_div (z, x, z, r);
  /* Note: if both x and y are 0, z is NAN, but the test below will
     be false. So, everything is fine. */
  if (mpfr_cmp_si (z, -1) < 0 || mpfr_cmp_ui (z, 1) > 0)
    {
      printf ("Error, -1 <= x/sqrt(x^2+y^2) <= 1 does not hold for r=%s\n",
              mpfr_print_rnd_mode (r));
      printf ("x="); mpfr_dump (x);
      printf ("y="); mpfr_dump (y);
      printf ("got "); mpfr_dump (z);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);
  mpfr_clear (t);
}

/* check sqrt(x^2) = x */
static void
test_property2 (mpfr_prec_t p, mpfr_rnd_t r)
{
  mpfr_t x, y;

  mpfr_init2 (x, p);
  mpfr_init2 (y, p);

  mpfr_urandomb (x, RANDS);
  mpfr_mul (y, x, x, r);
  mpfr_sqrt (y, y, r);
  if (mpfr_cmp (y, x))
    {
      printf ("Error, sqrt(x^2) = x does not hold for r=%s\n",
              mpfr_print_rnd_mode (r));
      printf ("x="); mpfr_dump (x);
      printf ("got "); mpfr_dump (y);
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

/* Bug reported by Fredrik Johansson, occurring when:
   - the precision of the result is a multiple of the number of bits
     per word (GMP_NUMB_BITS),
   - the rounding mode is to nearest (MPFR_RNDN),
   - internally, the result has to be rounded up to a power of 2.
*/
static void
bug20160120 (void)
{
  mpfr_t x, y;

  mpfr_init2 (x, 4 * GMP_NUMB_BITS);
  mpfr_init2 (y, GMP_NUMB_BITS);

  mpfr_set_ui (x, 1, MPFR_RNDN);
  mpfr_nextbelow (x);
  mpfr_sqrt (y, x, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_check (y));
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0);

  mpfr_set_prec (y, 2 * GMP_NUMB_BITS);
  mpfr_sqrt (y, x, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_check (y));
  MPFR_ASSERTN(mpfr_cmp_ui (y, 1) == 0);

  mpfr_clear(x);
  mpfr_clear(y);
}

/* Bug in mpfr_sqrt2 when prec(u) = 2*GMP_NUMB_BITS and the exponent of u is
   odd: the last bit of u is lost. */
static void
bug20160908 (void)
{
  mpfr_t r, u;
  int n = GMP_NUMB_BITS, ret;

  mpfr_init2 (r, 2*n - 1);
  mpfr_init2 (u, 2 * n);
  mpfr_set_ui_2exp (u, 1, 2*n-2, MPFR_RNDN); /* u=2^(2n-2) with exp(u)=2n-1 */
  mpfr_nextabove (u);
  /* now u = 2^(2n-2) + 1/2 */
  ret = mpfr_sqrt (r, u, MPFR_RNDZ);
  MPFR_ASSERTN(ret == -1 && mpfr_cmp_ui_2exp (r, 1, n-1) == 0);
  mpfr_clear (r);
  mpfr_clear (u);
}

static void
testall_rndf (mpfr_prec_t pmax)
{
  mpfr_t a, b, d;
  mpfr_prec_t pa, pb;

  for (pa = MPFR_PREC_MIN; pa <= pmax; pa++)
    {
      mpfr_init2 (a, pa);
      mpfr_init2 (d, pa);
      for (pb = MPFR_PREC_MIN; pb <= pmax; pb++)
        {
          mpfr_init2 (b, pb);
          mpfr_set_ui (b, 1, MPFR_RNDN);
          while (mpfr_cmp_ui (b, 4) < 0)
            {
              mpfr_sqrt (a, b, MPFR_RNDF);
              mpfr_sqrt (d, b, MPFR_RNDD);
              if (!mpfr_equal_p (a, d))
                {
                  mpfr_sqrt (d, b, MPFR_RNDU);
                  if (!mpfr_equal_p (a, d))
                    {
                      printf ("Error: mpfr_sqrt(a,b,RNDF) does not "
                              "match RNDD/RNDU\n");
                      printf ("b="); mpfr_dump (b);
                      printf ("a="); mpfr_dump (a);
                      exit (1);
                    }
                }
              mpfr_nextabove (b);
            }
          mpfr_clear (b);
        }
      mpfr_clear (a);
      mpfr_clear (d);
    }
}

/* test the case prec = GMP_NUMB_BITS */
static void
test_sqrt1n (void)
{
  mpfr_t r, u;
  int inex;

  MPFR_ASSERTD(GMP_NUMB_BITS >= 8); /* so that 15^2 is exactly representable */

  mpfr_init2 (r, GMP_NUMB_BITS);
  mpfr_init2 (u, GMP_NUMB_BITS);

  inex = mpfr_set_ui_2exp (u, 9 * 9, 2 * GMP_NUMB_BITS - 10, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  inex = mpfr_sqrt (r, u, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 9, GMP_NUMB_BITS - 5) == 0);

  inex = mpfr_set_ui_2exp (u, 15 * 15, 2 * GMP_NUMB_BITS - 10, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  inex = mpfr_sqrt (r, u, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 15, GMP_NUMB_BITS - 5) == 0);

  inex = mpfr_set_ui_2exp (u, 1, GMP_NUMB_BITS - 2, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  inex = mpfr_add_ui (u, u, 1, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  inex = mpfr_mul_2ui (u, u, GMP_NUMB_BITS, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  /* u = 2^(2*GMP_NUMB_BITS-2) + 2^GMP_NUMB_BITS, thus
     u = r^2 + 2^GMP_NUMB_BITS with r = 2^(GMP_NUMB_BITS-1).
     Should round to r+1 to nearest. */
  inex = mpfr_sqrt (r, u, MPFR_RNDN);
  MPFR_ASSERTN(inex > 0);
  inex = mpfr_sub_ui (r, r, 1, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, GMP_NUMB_BITS - 1) == 0);

  mpfr_clear (r);
  mpfr_clear (u);
}

static void
check_overflow (void)
{
  mpfr_t r, u;
  mpfr_prec_t p;
  mpfr_exp_t emax;
  int inex;

  emax = mpfr_get_emax ();
  for (p = MPFR_PREC_MIN; p <= 1024; p++)
    {
      mpfr_init2 (r, p);
      mpfr_init2 (u, p);

      set_emax (-1);
      mpfr_set_ui_2exp (u, 1, mpfr_get_emax () - 1, MPFR_RNDN);
      mpfr_nextbelow (u);
      mpfr_mul_2ui (u, u, 1, MPFR_RNDN);
      /* now u = (1 - 2^(-p))*2^emax is the largest number < +Inf,
         it square root is near 0.707 and has exponent 0 > emax */
      /* for RNDN, the result should be +Inf */
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_inf_p (r) && mpfr_sgn (r) > 0);
      /* for RNDA, the result should be +Inf */
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_inf_p (r) && mpfr_sgn (r) > 0);
      /* for RNDZ, the result should be u */
      inex = mpfr_sqrt (r, u, MPFR_RNDZ);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_equal_p (r, u));

      set_emax (0);
      mpfr_set_ui_2exp (u, 1, mpfr_get_emax () - 1, MPFR_RNDN);
      mpfr_nextbelow (u);
      mpfr_mul_2ui (u, u, 1, MPFR_RNDN);
      /* u = 1-2^(-p), its square root is > u, and should thus give +Inf when
         rounding away */
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_inf_p (r) && mpfr_sgn (r) > 0);

      mpfr_clear (r);
      mpfr_clear (u);
    }
  set_emax (emax);
}

static void
check_underflow (void)
{
  mpfr_t r, u;
  mpfr_prec_t p;
  mpfr_exp_t emin;
  int inex;

  emin = mpfr_get_emin ();
  for (p = MPFR_PREC_MIN; p <= 1024; p++)
    {
      mpfr_init2 (r, p);
      mpfr_init2 (u, p);

      set_emin (2);
      mpfr_set_ui_2exp (u, 1, mpfr_get_emin () - 1, MPFR_RNDN); /* u = 2 */
      /* for RNDN, since sqrt(2) is closer from 2 than 0, the result is 2 */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_equal_p (r, u));
      MPFR_ASSERTN(mpfr_underflow_p ());
      /* for RNDA, the result should be u, and there is underflow for p > 1,
         since for p=1 we have 1 < sqrt(2) < 2, but for p >= 2, sqrt(2) should
         be rounded to a number <= 1.5, which is representable */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_equal_p (r, u));
      MPFR_ASSERTN((p == 1 && !mpfr_underflow_p ()) ||
                   (p != 1 && mpfr_underflow_p ()));
      /* for RNDZ, the result should be +0 */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDZ);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_zero_p (r) && mpfr_signbit (r) == 0);
      MPFR_ASSERTN(mpfr_underflow_p ());

      /* generate an input u such that sqrt(u) < 0.5*2^emin but there is no
         underflow since sqrt(u) >= pred(0.5*2^emin), thus u >= 2^(2emin-2) */
      mpfr_set_ui_2exp (u, 1, 2 * mpfr_get_emin () - 2, MPFR_RNDN);
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex == 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex == 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDZ);
      MPFR_ASSERTN(inex == 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());

      /* next number */
      mpfr_set_ui_2exp (u, 1, 2 * mpfr_get_emin () - 2, MPFR_RNDN);
      mpfr_nextabove (u);
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex > 0);
      mpfr_nextbelow (r);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDZ);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());

      /* previous number */
      mpfr_set_ui_2exp (u, 1, 2 * mpfr_get_emin () - 2, MPFR_RNDN);
      mpfr_nextbelow (u);
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      /* since sqrt(u) is just below the middle between 0.5*2^emin and
         the previous number (with unbounded exponent range), there is
         underflow */
      MPFR_ASSERTN(mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDA);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(!mpfr_underflow_p ());
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDZ);
      MPFR_ASSERTN(inex < 0);
      mpfr_nextabove (r);
      MPFR_ASSERTN(mpfr_cmp_ui_2exp (r, 1, mpfr_get_emin () - 1) == 0);
      MPFR_ASSERTN(mpfr_underflow_p ());

      set_emin (3);
      mpfr_set_ui_2exp (u, 1, mpfr_get_emin () - 1, MPFR_RNDN); /* u = 4 */
      /* sqrt(u) = 2 = 0.5^2^(emin-1) should be rounded to +0 */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_zero_p (r) && mpfr_signbit (r) == 0);
      MPFR_ASSERTN(mpfr_underflow_p ());

      /* next number */
      mpfr_set_ui_2exp (u, 1, mpfr_get_emin () - 1, MPFR_RNDN); /* u = 4 */
      mpfr_nextabove (u);
      /* sqrt(u) should be rounded to 4 */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex > 0);
      MPFR_ASSERTN(mpfr_cmp_ui (r, 4) == 0);
      MPFR_ASSERTN(mpfr_underflow_p ());

      set_emin (4);
      mpfr_set_ui_2exp (u, 1, mpfr_get_emin () - 1, MPFR_RNDN); /* u = 8 */
      /* sqrt(u) should be rounded to +0 */
      mpfr_clear_flags ();
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      MPFR_ASSERTN(inex < 0);
      MPFR_ASSERTN(mpfr_zero_p (r) && mpfr_signbit (r) == 0);
      MPFR_ASSERTN(mpfr_underflow_p ());

      mpfr_clear (r);
      mpfr_clear (u);
    }
  set_emin (emin);
}

static void
coverage (void)
{
  mpfr_t r, t, u, v, w;
  mpfr_prec_t p;
  int inex;

  /* exercise even rule */
  for (p = MPFR_PREC_MIN; p <= 1024; p++)
    {
      mpfr_init2 (r, p);
      mpfr_init2 (t, p + 1);
      mpfr_init2 (u, 2 * p + 2);
      mpfr_init2 (v, p);
      mpfr_init2 (w, p);
      do
        mpfr_urandomb (v, RANDS);
      while (mpfr_zero_p (v));
      mpfr_set (w, v, MPFR_RNDN);
      mpfr_nextabove (w); /* w = nextabove(v) */
      mpfr_set (t, v, MPFR_RNDN);
      mpfr_nextabove (t);
      mpfr_sqr (u, t, MPFR_RNDN);
      inex = mpfr_sqrt (r, u, MPFR_RNDN);
      if (mpfr_min_prec (v) < p) /* v is even */
        {
          MPFR_ASSERTN(inex < 0);
          MPFR_ASSERTN(mpfr_equal_p (r, v));
        }
      else /* v is odd */
        {
          MPFR_ASSERTN(inex > 0);
          MPFR_ASSERTN(mpfr_equal_p (r, w));
        }
      mpfr_clear (r);
      mpfr_clear (t);
      mpfr_clear (u);
      mpfr_clear (v);
      mpfr_clear (w);
    }
}

#define TEST_FUNCTION test_sqrt
#define TEST_RANDOM_POS 8
#include "tgeneric.c"

int
main (void)
{
  mpfr_prec_t p;
  int k;

  tests_start_mpfr ();

  coverage ();
  check_underflow ();
  check_overflow ();
  testall_rndf (16);
  for (p = MPFR_PREC_MIN; p <= 128; p++)
    {
      test_property1 (p, MPFR_RNDN, 0);
      test_property1 (p, MPFR_RNDU, 0);
      test_property1 (p, MPFR_RNDA, 0);
      test_property1 (p, MPFR_RNDN, 1);
      test_property1 (p, MPFR_RNDU, 1);
      test_property1 (p, MPFR_RNDA, 1);
      test_property2 (p, MPFR_RNDN);
    }

  check_diverse ("635030154261163106768013773815762607450069292760790610550915652722277604820131530404842415587328", 160, "796887792767063979679855997149887366668464780637");
  special ();
  check_singular ();

  for (p=2; p<200; p++)
    for (k=0; k<200; k++)
      check_inexact (p);
  check_float();

  check3 ("-0.0", MPFR_RNDN, "0.0");
  check4 ("6.37983013646045901440e+32", MPFR_RNDN, "5.9bc5036d09e0c@13");
  check4 ("1.0", MPFR_RNDN, "1");
  check4 ("1.0", MPFR_RNDZ, "1");
  check4 ("3.725290298461914062500000e-9", MPFR_RNDN, "4@-4");
  check4 ("3.725290298461914062500000e-9", MPFR_RNDZ, "4@-4");
  check4 ("1190456976439861.0", MPFR_RNDZ, "2.0e7957873529a@6");
  check4 ("1219027943874417664.0", MPFR_RNDZ, "4.1cf2af0e6a534@7");
  /* the following examples are bugs in Cygnus compiler/system, found by
     Fabrice Rouillier while porting mpfr to Windows */
  check4 ("9.89438396044940256501e-134", MPFR_RNDU, "8.7af7bf0ebbee@-56");
  check4 ("7.86528588050363751914e+31", MPFR_RNDZ, "1.f81fc40f32062@13");
  check4 ("0.99999999999999988897", MPFR_RNDN, "f.ffffffffffff8@-1");
  check4 ("1.00000000000000022204", MPFR_RNDN, "1");
  /* the following examples come from the paper "Number-theoretic Test
   Generation for Directed Rounding" from Michael Parks, Table 4 */

  check4 ("78652858805036375191418371571712.0", MPFR_RNDN,
          "1.f81fc40f32063@13");
  check4 ("38510074998589467860312736661504.0", MPFR_RNDN,
          "1.60c012a92fc65@13");
  check4 ("35318779685413012908190921129984.0", MPFR_RNDN,
          "1.51d17526c7161@13");
  check4 ("26729022595358440976973142425600.0", MPFR_RNDN,
          "1.25e19302f7e51@13");
  check4 ("22696567866564242819241453027328.0", MPFR_RNDN,
          "1.0ecea7dd2ec3d@13");
  check4 ("22696888073761729132924856434688.0", MPFR_RNDN,
          "1.0ecf250e8e921@13");
  check4 ("36055652513981905145251657416704.0", MPFR_RNDN,
          "1.5552f3eedcf33@13");
  check4 ("30189856268896404997497182748672.0", MPFR_RNDN,
          "1.3853ee10c9c99@13");
  check4 ("36075288240584711210898775080960.0", MPFR_RNDN,
          "1.556abe212b56f@13");
  check4 ("72154663483843080704304789585920.0", MPFR_RNDN,
          "1.e2d9a51977e6e@13");

  check4 ("78652858805036375191418371571712.0", MPFR_RNDZ,
          "1.f81fc40f32062@13");
  check4 ("38510074998589467860312736661504.0", MPFR_RNDZ,
          "1.60c012a92fc64@13");
  check4 ("35318779685413012908190921129984.0", MPFR_RNDZ, "1.51d17526c716@13");
  check4 ("26729022595358440976973142425600.0", MPFR_RNDZ, "1.25e19302f7e5@13");
  check4 ("22696567866564242819241453027328.0", MPFR_RNDZ,
          "1.0ecea7dd2ec3c@13");
  check4 ("22696888073761729132924856434688.0", MPFR_RNDZ, "1.0ecf250e8e92@13");
  check4 ("36055652513981905145251657416704.0", MPFR_RNDZ,
          "1.5552f3eedcf32@13");
  check4 ("30189856268896404997497182748672.0", MPFR_RNDZ,
          "1.3853ee10c9c98@13");
  check4 ("36075288240584711210898775080960.0", MPFR_RNDZ,
          "1.556abe212b56e@13");
  check4 ("72154663483843080704304789585920.0", MPFR_RNDZ,
          "1.e2d9a51977e6d@13");

  check4 ("78652858805036375191418371571712.0", MPFR_RNDU,
          "1.f81fc40f32063@13");
  check4 ("38510074998589467860312736661504.0", MPFR_RNDU,
          "1.60c012a92fc65@13");
  check4 ("35318779685413012908190921129984.0", MPFR_RNDU,
          "1.51d17526c7161@13");
  check4 ("26729022595358440976973142425600.0", MPFR_RNDU,
          "1.25e19302f7e51@13");
  check4 ("22696567866564242819241453027328.0", MPFR_RNDU,
          "1.0ecea7dd2ec3d@13");
  check4 ("22696888073761729132924856434688.0", MPFR_RNDU,
          "1.0ecf250e8e921@13");
  check4 ("36055652513981905145251657416704.0", MPFR_RNDU,
          "1.5552f3eedcf33@13");
  check4 ("30189856268896404997497182748672.0", MPFR_RNDU,
          "1.3853ee10c9c99@13");
  check4 ("36075288240584711210898775080960.0", MPFR_RNDU,
          "1.556abe212b56f@13");
  check4 ("72154663483843080704304789585920.0", MPFR_RNDU,
          "1.e2d9a51977e6e@13");

  check4 ("78652858805036375191418371571712.0", MPFR_RNDD,
          "1.f81fc40f32062@13");
  check4 ("38510074998589467860312736661504.0", MPFR_RNDD,
          "1.60c012a92fc64@13");
  check4 ("35318779685413012908190921129984.0", MPFR_RNDD, "1.51d17526c716@13");
  check4 ("26729022595358440976973142425600.0", MPFR_RNDD, "1.25e19302f7e5@13");
  check4 ("22696567866564242819241453027328.0", MPFR_RNDD,
          "1.0ecea7dd2ec3c@13");
  check4 ("22696888073761729132924856434688.0", MPFR_RNDD, "1.0ecf250e8e92@13");
  check4 ("36055652513981905145251657416704.0", MPFR_RNDD,
          "1.5552f3eedcf32@13");
  check4 ("30189856268896404997497182748672.0", MPFR_RNDD,
          "1.3853ee10c9c98@13");
  check4 ("36075288240584711210898775080960.0", MPFR_RNDD,
          "1.556abe212b56e@13");
  check4 ("72154663483843080704304789585920.0", MPFR_RNDD,
          "1.e2d9a51977e6d@13");

  /* check that rounding away is just rounding toward positive infinity */
  check4 ("72154663483843080704304789585920.0", MPFR_RNDA,
          "1.e2d9a51977e6e@13");

  test_generic (MPFR_PREC_MIN, 300, 15);
  data_check ("data/sqrt", mpfr_sqrt, "mpfr_sqrt");
  bad_cases (mpfr_sqrt, mpfr_sqr, "mpfr_sqrt", 0, -256, 255, 4, 128, 80, 50);

  bug20160120 ();
  bug20160908 ();
  test_sqrt1n ();

  tests_end_mpfr ();
  return 0;
}
