/* Test file for mpfr_subnormalize.

Copyright 2005-2025 Free Software Foundation, Inc.
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

static const struct {
  const char *in;
  int i;
  mpfr_rnd_t rnd;
  const char *out;
  int j;
} tab[] = {
  /* 2nd field: old_inex value */
  /* 4th field: use the mpfr_dump format, in case of error. */
  /* 5th field: inexact field */
  {"1E1",  0, MPFR_RNDN, "0.100000000E2", 0},
  {"1E1", -1, MPFR_RNDZ, "0.100000000E2", -1},
  {"1E1", -1, MPFR_RNDD, "0.100000000E2", -1},
  {"1E1",  1, MPFR_RNDU, "0.100000000E2", 1},
  {"0.10000E-10", 0, MPFR_RNDN, "0.100000000E-10", 0},
  {"0.10001E-10", 0, MPFR_RNDN, "0.100000000E-10", -1},
  {"0.11001E-10", 0, MPFR_RNDN, "0.100000000E-9", 1},
  {"0.11001E-10", 0, MPFR_RNDZ, "0.100000000E-10", -1},
  {"0.11001E-10", 0, MPFR_RNDU, "0.100000000E-9", 1},
  {"0.11000E-10", 0, MPFR_RNDN, "0.100000000E-9", 1},
  {"0.11000E-10", -1, MPFR_RNDN, "0.100000000E-9", 1},
  {"0.11000E-10", 1, MPFR_RNDN, "0.100000000E-10", -1},
  {"0.11111E-8", 0, MPFR_RNDN, "0.100000000E-7", 1},
  {"0.10111E-8", 0, MPFR_RNDN, "0.110000000E-8", 1},
  {"0.11110E-8", -1, MPFR_RNDN, "0.100000000E-7", 1},
  {"0.10110E-8", 1, MPFR_RNDN, "0.101000000E-8", -1}
};

static void
check1 (void)
{
  mpfr_t x;
  int i, j, k, s, old_inex, tiny, expj;
  mpfr_exp_t emin, emax;
  unsigned int expflags, flags;

  emin = mpfr_get_emin ();
  emax = mpfr_get_emax ();

  mpfr_set_default_prec (9);
  set_emin (-10);
  set_emax (10);

  mpfr_init (x);
  for (i = 0; i < numberof (tab); i++)
    for (s = 0; s <= (tab[i].rnd == MPFR_RNDN); s++)
      for (k = 0; k <= 1; k++)
        {
          mpfr_set_str (x, tab[i].in, 2, MPFR_RNDN);
          old_inex = tab[i].i;
          expj = tab[i].j;
          if (s)
            {
              mpfr_neg (x, x, MPFR_RNDN);
              old_inex = - old_inex;
              expj = - expj;
            }
          if (k && old_inex)
            old_inex = old_inex < 0 ? INT_MIN : INT_MAX;
          tiny = MPFR_GET_EXP (x) <= -3;
          mpfr_clear_flags ();
          j = mpfr_subnormalize (x, old_inex, tab[i].rnd);
          expflags =
            (tiny ? MPFR_FLAGS_UNDERFLOW : 0) |
            (expj ? MPFR_FLAGS_INEXACT : 0);
          flags = __gmpfr_flags;
          if (s)
            mpfr_neg (x, x, MPFR_RNDN);
          if (mpfr_cmp_str (x, tab[i].out, 2, MPFR_RNDN) != 0 ||
              flags != expflags || ! SAME_SIGN (j, expj))
            {
              const char *sgn = s ? "-" : "";
              printf ("Error for i = %d (old_inex = %d), k = %d, x = %s%s\n"
                      "Expected: %s%s\nGot:      ", i, old_inex, k,
                      sgn, tab[i].in, sgn, tab[i].out);
              if (s)
                mpfr_neg (x, x, MPFR_RNDN);
              mpfr_dump (x);
              printf ("Expected flags = %u, got %u\n", expflags, flags);
              printf ("Expected ternary value = %d, got %d\n", expj, j);
              exit (1);
            }
        }
  mpfr_clear (x);

  MPFR_ASSERTN (mpfr_get_emin () == -10);
  MPFR_ASSERTN (mpfr_get_emax () == 10);

  set_emin (emin);
  set_emax (emax);
}

/* bug found by Kevin P. Rauch on 22 Oct 2007 */
static void
check2 (void)
{
  mpfr_t x, y, z;
  int tern;
  mpfr_exp_t emin;

  emin = mpfr_get_emin ();

  mpfr_init2 (x, 32);
  mpfr_init2 (y, 32);
  mpfr_init2 (z, 32);

  mpfr_set_ui (x, 0xC0000000U, MPFR_RNDN);
  mpfr_neg (x, x, MPFR_RNDN);
  mpfr_set_ui (y, 0xFFFFFFFEU, MPFR_RNDN);
  mpfr_set_exp (x, 0);
  mpfr_set_exp (y, 0);
  set_emin (-29);

  tern = mpfr_mul (z, x, y, MPFR_RNDN);
  /* z = -0.BFFFFFFE, tern > 0 */

  tern = mpfr_subnormalize (z, tern, MPFR_RNDN);
  /* z should be -0.75 */
  MPFR_ASSERTN (tern < 0 && mpfr_cmp_si_2exp (z, -3, -2) == 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);

  MPFR_ASSERTN (mpfr_get_emin () == -29);

  set_emin (emin);
}

/* bug found by Kevin P. Rauch on 22 Oct 2007 */
static void
check3 (void)
{
  mpfr_t x, y, z;
  int tern;
  mpfr_exp_t emin;

  emin = mpfr_get_emin ();

  mpfr_init2 (x, 32);
  mpfr_init2 (y, 32);
  mpfr_init2 (z, 32);

  mpfr_set_ui (x, 0xBFFFFFFFU, MPFR_RNDN); /* 3221225471/2^32 */
  mpfr_set_ui (y, 0x80000001U, MPFR_RNDN); /* 2147483649/2^32 */
  mpfr_set_exp (x, 0);
  mpfr_set_exp (y, 0);
  set_emin (-1);

  /* the exact product is 6917529028714823679/2^64, which is rounded to
     3/8 = 0.375, which is smaller, thus tern < 0 */
  tern = mpfr_mul (z, x, y, MPFR_RNDN);
  MPFR_ASSERTN (tern < 0 && mpfr_cmp_ui_2exp (z, 3, -3) == 0);

  tern = mpfr_subnormalize (z, tern, MPFR_RNDN);
  /* since emin = -1, and EXP(z)=-1, z should be rounded to precision
     EXP(z)-emin+1 = 1, i.e., z should be a multiple of the smallest possible
     positive representable value with emin=-1, which is 1/4. The two
     possible values are 1/4 and 2/4, which are at equal distance of z.
     But since tern < 0, we should choose the largest value, i.e., 2/4. */
  MPFR_ASSERTN (tern > 0 && mpfr_cmp_ui_2exp (z, 1, -1) == 0);

  /* here is another test for the alternate case, where z was rounded up
     first, thus we have to round down */
  mpfr_set_str_binary (x, "0.11111111111010110101011011011011");
  mpfr_set_str_binary (y, "0.01100000000001111100000000001110");
  tern = mpfr_mul (z, x, y, MPFR_RNDN);
  MPFR_ASSERTN (tern > 0 && mpfr_cmp_ui_2exp (z, 3, -3) == 0);
  tern = mpfr_subnormalize (z, tern, MPFR_RNDN);
  MPFR_ASSERTN (tern < 0 && mpfr_cmp_ui_2exp (z, 1, -2) == 0);

  /* finally the case where z was exact, which we simulate here */
  mpfr_set_ui_2exp (z, 3, -3, MPFR_RNDN);
  tern = mpfr_subnormalize (z, 0, MPFR_RNDN);
  MPFR_ASSERTN (tern > 0 && mpfr_cmp_ui_2exp (z, 1, -1) == 0);

  mpfr_clear (x);
  mpfr_clear (y);
  mpfr_clear (z);

  MPFR_ASSERTN (mpfr_get_emin () == -1);

  set_emin (emin);
}

/* check mpfr_subnormize with RNDNA (experimental) */
static void
check_rndna (void)
{
  mpfr_t x, y;
  int inex;
  mpfr_exp_t emin = mpfr_get_emin ();

  mpfr_init2 (x, 11);
  mpfr_init2 (y, 9);

  mpfr_set_str_binary (x, "0.1111111010E-14");
  inex = mpfr_set (y, x, MPFR_RNDN);
  MPFR_ASSERTN(inex == 0);
  set_emin (-21);
  inex = mpfr_subnormalize (y, inex, MPFR_RNDNA);
  /* mpfr_subnormalize rounds to precision EXP(y)-emin+1 = 8,
     thus should round to 0.111111110E-14 */
  mpfr_set_str_binary (x, "0.111111110E-14");
  MPFR_ASSERTN(mpfr_cmp (y, x) == 0);
  MPFR_ASSERTN(inex > 0);

  /* now consider x slightly larger: we should get the same result */
  mpfr_set_str_binary (x, "0.1111111011E-14");
  inex = mpfr_set (y, x, MPFR_RNDN);
  MPFR_ASSERTN(inex > 0);
  inex = mpfr_subnormalize (y, inex, MPFR_RNDNA);
  mpfr_set_str_binary (x, "0.111111110E-14");
  MPFR_ASSERTN(mpfr_cmp (y, x) == 0);
  MPFR_ASSERTN(inex > 0);

  /* now consider x slightly smaller: we should get a different result */
  mpfr_set_str_binary (x, "0.11111110001E-14");
  inex = mpfr_set (y, x, MPFR_RNDN);
  MPFR_ASSERTN(inex < 0);
  inex = mpfr_subnormalize (y, inex, MPFR_RNDNA);
  mpfr_set_str_binary (x, "0.111111100E-14");
  MPFR_ASSERTN(mpfr_cmp (y, x) == 0);
  MPFR_ASSERTN(inex < 0);

  mpfr_clear (x);
  mpfr_clear (y);
  set_emin (emin);
}

/* exercise a corner case of mpfr_subnormalize:
   y = 1xxx...xxx0|100000| with old_inexact = -1 */
static void
coverage (void)
{
  mpfr_t y;
  int inex;

  mpfr_init2 (y, 42);
  mpfr_set_ui_2exp (y, 131073, mpfr_get_emin () - 2, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_get_exp (y) == mpfr_get_emin () + 16);
  /* mpfr_subnormalize rounds y to precision EXP(y) - emin + 1, thus 17 */
  inex = mpfr_subnormalize (y, -1, MPFR_RNDN);
  MPFR_ASSERTN(mpfr_cmp_ui_2exp (y, 65537, mpfr_get_emin () - 1) == 0);
  MPFR_ASSERTN(inex > 0);
  mpfr_clear (y);
}

/* check underflow is after rounding */
static void
check_underflow (void)
{
  mpfr_t zm;
  mpfr_exp_t emin = mpfr_get_emin ();
  int i;
  int Inex[] = { -1, 0, +1 };

  set_emin (-1073);
  mpfr_init2 (zm, 53);
  for (i = 0; i < 3; i++)
    {
      mpfr_flags_clear (MPFR_FLAGS_UNDERFLOW);
      mpfr_set_str (zm, "0x3.ffffffffffffep-1024", 16, MPFR_RNDN);
      mpfr_subnormalize (zm, Inex[i], MPFR_RNDU);
      if (!mpfr_flags_test (MPFR_FLAGS_UNDERFLOW))
        {
          printf ("Missing underflow for zm=0x3.ffffffffffffep-1024\n");
          printf ("old_inex = %d\n", Inex[i]);
          printf ("after mpfr_subnormalize, zm=");
          mpfr_dump (zm);
          printf ("\n");
          exit (1);
        }

      mpfr_flags_clear (MPFR_FLAGS_UNDERFLOW);
      mpfr_set_str (zm, "0x1p-1022", 16, MPFR_RNDN);
      mpfr_subnormalize (zm, Inex[i], MPFR_RNDU);
      if (mpfr_flags_test (MPFR_FLAGS_UNDERFLOW))
        {
          printf ("Spurious underflow for zm=0x1p-1022\n");
          printf ("old_inex = %d\n", Inex[i]);
          printf ("after mpfr_subnormalize, zm=");
          mpfr_dump (zm);
          printf ("\n");
          exit (1);
        }
    }

  mpfr_clear (zm);
  set_emin (emin);
}

int
main (int argc, char *argv[])
{
  tests_start_mpfr ();

  coverage ();
  check_rndna ();
  check1 ();
  check2 ();
  check3 ();
  check_underflow ();

  tests_end_mpfr ();
  return 0;
}
