/* Test file for mpfr_nrandom

Copyright 2011-2026 Free Software Foundation, Inc.
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

/* The number of variants of nrandom */
#define NRANDOM_VERSIONS 2

static void
test_special (int version, mpfr_prec_t p)
{
  mpfr_t x;
  int inexact;

  mpfr_init2 (x, p);

  inexact = (version == 1 ? mpfr_nrandom_v1 : mpfr_nrandom_v2)
    (x, RANDS, MPFR_RNDN);
  if (inexact == 0)
    {
      printf ("Error: mpfr_nrandom_v%d() returns a zero ternary value.\n",
              version);
      exit (1);
    }

  mpfr_clear (x);
}

#define NRES 10

/* First NRES entries for v1, second NRES entries for v2 */
static const char *res[2*NRES] = {
  "-2.07609e2d96da78b2d6bea3ab30d4359222a82f8a35e4e4464303ad4808f57458@0",
  "1.a4650a963ab9f266ed009ee96c8788f6b88212f5f2a4d4aef65db2a9e57c44bc@-1",
  "c.0b3cda7f370a36febed972dbb47f2503f7e08a651edbf12d0303d968257841b0@-1",
  "-7.4dffd19868c4bc2ec5b6c8311e0f190e1c97b575b7fdabec897e5de2a1d802b8@-1",
  "1.40f6204ded71a4346ed17094863347b8c735e62712cc0b4d0c5402ee310d9714@0",
  "-3.09fd2a1fc234e23bfa048dbf1e7850ac6cbea2514a0f3ce011e964d9f331cbcc@-1",
  "-1.104e769aadb5fce5a7ad1c546e91b889829a76920c7cc7ac4cbd12009451ce90@0",
  "3.0a08181e342b02187463c0025f895b41ddb7076c5bf157e3b898e9248baf4ad4@-1",
  "-d.44fda7a51276b722ebc88dd016b7d9d7ea5ba682282a42cdef6948312e5dcf70@-1",
  "1.5bf69aff31bb3e6430cc263fdd45ef2c70a779984e764524bc35a9cb4a430dd0@0",

  "8.90acc8b8a980ab5564006cb05bf8573b7bc3ed6a06431d4ce094858d11fa8dd0@-1",
  "-1.d025e3f430d4359222a82f8a35e4e4464303ad4808f5745781361c0dd614fc2e@0",
  "-e.a66fb696b88212f5f2a4d4aef65db2a9e57c44bc58e39520eda7e001a4650aa0@-1",
  "1.80f0a9628257841b744575cc65865e54f4cf43f4f1527ebd09f395470195253e@0",
  "5.b7fdabe55c5dd92d768983f5fe84a7c8be40c10868c4bc2ec5b6c8311e0f1910@-1",
  "a.1173efbbb17b798626f226dde71006d0de654a2523935c07da17dfd40f6204d0@-1",
  "-8.ed5bb0dab7e39130db211bc4bbce42f0fd4f40db3bc6be9e87c779358e4afc40@-1",
  "3.09fd2a1bfa048dbf1e7850ac6cbea2514a0f3ce011e964d9f331cbcab2efdf10@-1",
  "1.21209642062a049a1ad0e1792dc8a971020cab00adb5fce5a7ad1c546e91b88a@0",
  "3.0a0818187463c0025f895b41ddb7076c5bf157e3b898e9248baf4ad266fc2f70@-1" };

/* If checkval is true, check the obtained results by using a fixed seed
   for reproducibility. */
static void
test_nrandom (int version, long nbtests, mpfr_prec_t prec, mpfr_rnd_t rnd,
              int verbose, int checkval)
{
  gmp_randstate_t s;
  mpfr_t *t;
  int i, inexact;

  if (checkval)
    {
      gmp_randinit_default (s);
      gmp_randseed_ui (s, 17);
      nbtests = NRES;
    }

  t = (mpfr_t *) tests_allocate (nbtests * sizeof (mpfr_t));

  for (i = 0; i < nbtests; ++i)
    mpfr_init2 (t[i], prec);

  for (i = 0; i < nbtests; i++)
    {
      inexact = (version == 1 ? mpfr_nrandom_v1 : mpfr_nrandom_v2)
        (t[i], checkval ? s : RANDS, MPFR_RNDN);

      if (checkval &&
          mpfr_cmp_str (t[i], res[NRES*(version-1) + i], 16, MPFR_RNDN) != 0)
        {
          printf ("Unexpected value in test_nrandom_v%d().\n"
                  "Expected %s\n"
                  "Got      ", version, res[NRES*(version-1) + i]);
          mpfr_out_str (stdout, 16, 0, t[i], MPFR_RNDN);
          printf ("\n");
          exit (1);
        }

      if (inexact == 0)
        {
          /* one call in the loop pretended to return an exact number! */
          printf ("Error: mpfr_nrandom_v%d() returns a zero ternary value.\n",
                  version);
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

      mpfr_printf ("Average v%d = %.5Rf\nVariance v%d = %.5Rf\n",
                   version, av, version, va);
      mpfr_clear (av);
      mpfr_clear (va);
      mpfr_clear (tmp);
    }
#endif /* HAVE_STDARG */

  for (i = 0; i < nbtests; ++i)
    mpfr_clear (t[i]);
  tests_free (t, nbtests * sizeof (mpfr_t));
  if (checkval)
    gmp_randclear (s);
  return;
}

int
main (int argc, char *argv[])
{
  long nbtests;
  int verbose, v;

  tests_start_mpfr ();

  verbose = 0;
  nbtests = 10;
  if (argc > 1)
    {
      /* Number of values in argument. Note that the mpfr_clear loop above
         is in O(n^2) until the FIXME for tests_memory_find() in memory.c
         is resolved (the search in tests_memory_find() is in O(n), while
         it could be in almost constant time). */
      long a = atol (argv[1]);
      verbose = 1;
      if (a != 0)
        nbtests = a;
    }

  for (v = 1; v <= NRANDOM_VERSIONS; ++v)
    {
      test_nrandom (v, nbtests, 420, MPFR_RNDN, verbose, 0);

#ifndef MPFR_USE_MINI_GMP
      /* The random generator in mini-gmp is not deterministic. */
      test_nrandom (v, 0, 256, MPFR_RNDN, 0, 1);
#endif /* MPFR_USE_MINI_GMP */

      test_special (v, 2);
      test_special (v, 42000);
    }

  tests_end_mpfr ();
  return 0;
}
