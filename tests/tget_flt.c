/* Test file for mpfr_get_flt and mpfr_set_flt

Copyright 2009-2025 Free Software Foundation, Inc.
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

#include <float.h>     /* for FLT_MIN */

#include "mpfr-test.h"
#include "ieee_floats.h"

/* return non-zero iff f == g, and if both are zero check the sign */
static int
equal_flt (float f, float g)
{
  if (f != 0.0)
    return f == g;
  else if (g != 0)
    return 0;
  else /* f = g = 0: check they have the same sign */
    {
      int sf, sg;
      mpfr_t z;

      mpfr_init2 (z, MPFR_PREC_MIN);
      mpfr_set_flt (z, f, MPFR_RNDN);
      sf = mpfr_signbit (z);
      mpfr_set_flt (z, g, MPFR_RNDN);
      sg = mpfr_signbit (z);
      mpfr_clear (z);
      return !sf == !sg;
    }
}

int
main (void)
{
  mpfr_t x, y;
  float f, g;
  int i;
#if !defined(MPFR_ERRDIVZERO)
  float infp;
#endif

  tests_start_mpfr ();

#if !defined(MPFR_ERRDIVZERO)
  infp = (float) MPFR_DBL_INFP;
  if (infp * 0.5 != infp)
    {
      fprintf (stderr, "Error, FLT_MAX + FLT_MAX does not yield INFP\n");
      fprintf (stderr, "(this is probably a compiler bug, please report)\n");
      exit (1);
    }
#endif

  mpfr_init2 (x, 24);
  mpfr_init2 (y, 24);

#if !defined(MPFR_ERRDIVZERO)
  mpfr_set_nan (x);
  f = mpfr_get_flt (x, MPFR_RNDN);
  if (! DOUBLE_ISNAN (f))
    {
      printf ("Error for mpfr_get_flt(NaN)\n");
      printf ("got f=%f\n", f);
      exit (1);
    }
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (! mpfr_nan_p (x))
    {
      printf ("Error for mpfr_set_flt(NaN)\n");
      exit (1);
    }

  mpfr_set_inf (x, 1);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (! mpfr_inf_p (x) || mpfr_sgn (x) < 0)
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(+Inf)):\n");
      printf ("f=%f, expected -Inf\n", f);
      printf ("got "); mpfr_dump (x);
      exit (1);
    }

  mpfr_set_inf (x, -1);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (! mpfr_inf_p (x) || mpfr_sgn (x) > 0)
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(-Inf)):\n");
      printf ("f=%f, expected -Inf\n", f);
      printf ("got "); mpfr_dump (x);
      exit (1);
    }
#endif

  mpfr_set_ui (x, 0, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (! mpfr_zero_p (x) || MPFR_IS_NEG (x))
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(+0))\n");
      exit (1);
    }

#ifdef HAVE_SIGNEDZ
  mpfr_set_ui (x, 0, MPFR_RNDN);
  mpfr_neg (x, x, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (! mpfr_zero_p (x) || MPFR_IS_POS (x))
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(-0))\n");
      exit (1);
    }
#endif  /* HAVE_SIGNEDZ */

  mpfr_set_ui (x, 17, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (mpfr_cmp_ui (x, 17) != 0)
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(17))\n");
      printf ("expected 17\n");
      printf ("got      ");
      mpfr_dump (x);
      exit (1);
    }

  mpfr_set_si (x, -42, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN);
  mpfr_set_flt (x, f, MPFR_RNDN);
  if (mpfr_cmp_si (x, -42) != 0)
    {
      printf ("Error for mpfr_set_flt(mpfr_get_flt(-42))\n");
      printf ("expected -42\n");
      printf ("got      ");
      mpfr_dump (x);
      exit (1);
    }

  mpfr_set_si_2exp (x, 1, -126, MPFR_RNDN);
  for (i = -126; i < 128; i++)
    {
      f = mpfr_get_flt (x, MPFR_RNDN);
      mpfr_set_flt (y, f, MPFR_RNDN);
      if (mpfr_cmp (x, y) != 0)
        {
          printf ("Error for mpfr_set_flt(mpfr_get_flt(x))\n");
          printf ("expected "); mpfr_dump (x);
          printf ("got      "); mpfr_dump (y);
          exit (1);
        }
      mpfr_mul_2ui (x, x, 1, MPFR_RNDN);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_si_2exp (x, 1, -126, MPFR_RNDN);
  for (i = -126; i < 128; i++)
    {
      mpfr_nextbelow (x);
      f = mpfr_get_flt (x, MPFR_RNDN);
      mpfr_nextabove (x);
      mpfr_set_flt (y, f, MPFR_RNDN);
      if (mpfr_cmp (x, y) != 0)
        {
          printf ("Error for mpfr_set_flt(mpfr_get_flt(x))\n");
          printf ("expected "); mpfr_dump (x);
          printf ("got      "); mpfr_dump (y);
          exit (1);
        }
      mpfr_mul_2ui (x, x, 1, MPFR_RNDN);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_si_2exp (x, 1, -126, MPFR_RNDN);
  for (i = -126; i < 128; i++)
    {
      mpfr_nextabove (x);
      f = mpfr_get_flt (x, MPFR_RNDN);
      mpfr_nextbelow (x);
      mpfr_set_flt (y, f, MPFR_RNDN);
      if (mpfr_cmp (x, y) != 0)
        {
          printf ("Error for mpfr_set_flt(mpfr_get_flt(x))\n");
          printf ("expected "); mpfr_dump (x);
          printf ("got      "); mpfr_dump (y);
          exit (1);
        }
      mpfr_mul_2ui (x, x, 1, MPFR_RNDN);
    }

  if (have_subnorm_flt ())
    for (i = 0; i < 2; i++)
      {
        /* We assume here that the format of float is binary32, a.k.a.
           IEEE single precision (Annex F of ISO C for IEEE 754 support),
           as this is the case on all machines nowadays. Then 2^(-150) is
           halfway between 0 and the smallest positive float 2^(-149). */
        mpfr_set_si_2exp (x, 1, -150, MPFR_RNDN);
        g = 0.0;
        if (i == 1)
          {
            mpfr_neg (x, x, MPFR_RNDN);
            g = -g;
          }
        f = mpfr_get_flt (x, MPFR_RNDN);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-150),RNDN)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDZ);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-150),RNDZ)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, (i == 0) ? MPFR_RNDD : MPFR_RNDU);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-150),%s)\n",
                    i == 0 ? "RNDD" : "RNDU");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        g = FLT_MIN * FLT_EPSILON;
        if (i == 1)
          g = -g;
        f = mpfr_get_flt (x, (i == 0) ? MPFR_RNDU : MPFR_RNDD);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-150),%s)\n",
                    i == 0 ? "RNDU" : "RNDD");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDA);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-150),RNDA)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }

        mpfr_set_si_2exp (x, 1, -151, MPFR_RNDN);
        g = 0.0;
        if (i == 1)
          {
            mpfr_neg (x, x, MPFR_RNDN);
            g = -g;
          }
        f = mpfr_get_flt (x, MPFR_RNDN);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-151),RNDN)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDZ);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-151),RNDZ)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, (i == 0) ? MPFR_RNDD : MPFR_RNDU);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-151),%s)\n",
                    i == 0 ? "RNDD" : "RNDU");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        g = FLT_MIN * FLT_EPSILON;
        if (i == 1)
          g = -g;
        f = mpfr_get_flt (x, (i == 0) ? MPFR_RNDU : MPFR_RNDD);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-151),%s)\n",
                    i == 0 ? "RNDU" : "RNDD");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDA);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-151),RNDA)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }

        mpfr_set_si_2exp (x, 1, -149, MPFR_RNDN);
        g = FLT_MIN * FLT_EPSILON;
        if (i == 1)
          {
            mpfr_neg (x, x, MPFR_RNDN);
            g = -g;
          }
        f = mpfr_get_flt (x, MPFR_RNDN);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-149),RNDN)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDZ);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-149),RNDZ)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDD);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-149),RNDD)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDU);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-149),RNDU)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
        f = mpfr_get_flt (x, MPFR_RNDA);
        if (!equal_flt (f, g))
          {
            printf ("Error for mpfr_get_flt(2^(-149),RNDA)\n");
            printf ("expected %.8e, got %.8e\n", g, f);
            exit (1);
          }
      }  /* for loop with tests on subnormals */

  mpfr_set_si_2exp (x, 1, 128, MPFR_RNDN);
  g = FLT_MAX;
  f = mpfr_get_flt (x, MPFR_RNDZ);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128,RNDZ)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDD);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128,RNDD)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }

  /* same for negative x */
  mpfr_set_si_2exp (x, -1, 128, MPFR_RNDN);
  g = -FLT_MAX;
  f = mpfr_get_flt (x, MPFR_RNDZ);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(-2^128,RNDZ)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDU);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(-2^128,RNDD)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }

#if !defined(MPFR_ERRDIVZERO)
  mpfr_set_si_2exp (x, 1, 128, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN); /* 2^128 rounds to itself with extended
                                      exponent range, we should get +Inf */
  g = infp;
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128,RNDN)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDU);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128,RNDU)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDA);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128,RNDA)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }

  /* same for negative x */
  mpfr_set_si_2exp (x, -1, 128, MPFR_RNDN);
  f = mpfr_get_flt (x, MPFR_RNDN); /* -2^128 rounds to itself with extended
                                      exponent range, we should get +Inf */
  g = -infp;
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(-2^128,RNDN)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDD);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(-2^128,RNDD)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDA);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(-2^128,RNDA)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
#endif

  /* corner case: take x with 25 bits just below 2^128 */
  mpfr_set_prec (x, 25);
  mpfr_set_si_2exp (x, 1, 128, MPFR_RNDN);
  mpfr_nextbelow (x);
  g = FLT_MAX;
  f = mpfr_get_flt (x, MPFR_RNDZ);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128*(1-2^(-25)),RNDZ)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDD);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128*(1-2^(-25)),RNDD)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
#if !defined(MPFR_ERRDIVZERO)
  f = mpfr_get_flt (x, MPFR_RNDN); /* first round to 2^128 (even rule),
                                      thus we should get +Inf */
  g = infp;
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128*(1-2^(-25)),RNDN)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDU);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128*(1-2^(-25)),RNDU)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
  f = mpfr_get_flt (x, MPFR_RNDA);
  if (f != g)
    {
      printf ("Error for mpfr_get_flt(2^128*(1-2^(-25)),RNDA)\n");
      printf ("expected %.8e, got %.8e\n", g, f);
      exit (1);
    }
#endif

  mpfr_clear (x);
  mpfr_clear (y);

  tests_end_mpfr ();
  return 0;
}
