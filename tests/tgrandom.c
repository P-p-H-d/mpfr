/* Test file for mpfr_grandom

Copyright 2011-2025 Free Software Foundation, Inc.
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

#define _MPFR_NO_DEPRECATED_GRANDOM
#include "mpfr-test.h"

static void
test_special (mpfr_prec_t p)
{
  mpfr_t x;
  int inexact;

  mpfr_init2 (x, p);

  inexact = mpfr_grandom (x, NULL, RANDS, MPFR_RNDN);
  if (((unsigned int) inexact & 3) == 0)
    {
      printf ("Error: mpfr_grandom() returns a zero ternary value.\n");
      exit (1);
    }
  if (((unsigned int) inexact & (3 << 2)) != 0)
    {
      printf ("Error: the second ternary value of mpfr_grandom(x, NULL, ...)"
              " must be 0.\n");
      exit (1);
    }

  mpfr_clear (x);
}


static void
test_grandom (long nbtests, mpfr_prec_t prec, mpfr_rnd_t rnd,
              int verbose)
{
  mpfr_t *t;
  int i, inexact;

  nbtests = (nbtests & 1) ? (nbtests + 1) : nbtests;
  t = (mpfr_t *) tests_allocate (nbtests * sizeof (mpfr_t));

  for (i = 0; i < nbtests; ++i)
    mpfr_init2 (t[i], prec);

  for (i = 0; i < nbtests; i += 2)
    {
      inexact = mpfr_grandom (t[i], t[i + 1], RANDS, MPFR_RNDN);
      if (((unsigned int) inexact & 3) == 0 ||
          ((unsigned int) inexact & (3 << 2)) == 0)
        {
          /* one call in the loop pretended to return an exact number! */
          printf ("Error: mpfr_grandom() returns a zero ternary value.\n");
          exit (1);
        }
    }

#if defined(HAVE_STDARG) && !defined(MPFR_USE_MINI_GMP)
  if (verbose)
    {
      mpfr_t av, va, tmp;

      mpfr_init2 (av, prec);
      mpfr_init2 (va, prec);
      mpfr_init2 (tmp, prec);

      mpfr_set_ui (av, 0, MPFR_RNDN);
      mpfr_set_ui (va, 0, MPFR_RNDN);
      for (i = 0; i < nbtests; ++i)
        {
          mpfr_add (av, av, t[i], MPFR_RNDN);
          mpfr_sqr (tmp, t[i], MPFR_RNDN);
          mpfr_add (va, va, tmp, MPFR_RNDN);
        }
      mpfr_div_ui (av, av, nbtests, MPFR_RNDN);
      mpfr_div_ui (va, va, nbtests, MPFR_RNDN);
      mpfr_sqr (tmp, av, MPFR_RNDN);
      mpfr_sub (va, va, av, MPFR_RNDN);

      mpfr_printf ("Average = %.5Rf\nVariance = %.5Rf\n", av, va);
      mpfr_clear (av);
      mpfr_clear (va);
      mpfr_clear (tmp);
    }
#endif /* HAVE_STDARG */

  for (i = 0; i < nbtests; ++i)
    mpfr_clear (t[i]);
  tests_free (t, nbtests * sizeof (mpfr_t));
  return;
}


int
main (int argc, char *argv[])
{
  long nbtests;
  int verbose;
  int i;

  tests_start_mpfr ();

  verbose = 0;
  nbtests = 10;
  if (argc > 1)
    {
      long a = atol (argv[1]);
      verbose = 1;
      if (a != 0)
        nbtests = a;
    }

  test_grandom (nbtests, 420, MPFR_RNDN, verbose);
  test_special (2);
  test_special (42000);
  /* the following should exercise the case "Extend by 32 bits" in grandom.c */
  for (i = 0; i < 10000; i++)
    test_special (1);

  tests_end_mpfr ();
  return 0;
}
