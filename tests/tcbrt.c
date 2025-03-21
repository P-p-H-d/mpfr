/* Test file for mpfr_cbrt.

Copyright 2002-2025 Free Software Foundation, Inc.
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

static void
special (void)
{
  mpfr_t x, y;
  int inex;

  mpfr_init (x);
  mpfr_init (y);

  /* cbrt(NaN) = NaN */
  mpfr_set_nan (x);
  inex = mpfr_cbrt (y, x, MPFR_RNDN);
  if (!mpfr_nan_p (y))
    {
      printf ("Error: cbrt(NaN) is not NaN\n");
      exit (1);
    }
  if (inex != 0)
    {
      printf ("Error: cbrt(NaN): incorrect ternary value %d\n", inex);
      exit (1);
    }

  /* cbrt(+Inf) = +Inf */
  mpfr_set_inf (x, 1);
  inex = mpfr_cbrt (y, x, MPFR_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) < 0)
    {
      printf ("Error: cbrt(+Inf) <> +Inf\n");
      exit (1);
    }
  if (inex != 0)
    {
      printf ("Error: cbrt(+Inf): incorrect ternary value %d\n", inex);
      exit (1);
    }

  /* cbrt(-Inf) =  -Inf */
  mpfr_set_inf (x, -1);
  inex = mpfr_cbrt (y, x, MPFR_RNDN);
  if (!mpfr_inf_p (y) || mpfr_sgn (y) > 0)
    {
      printf ("Error: cbrt(-Inf) <> -Inf\n");
      exit (1);
    }
  if (inex != 0)
    {
      printf ("Error: cbrt(-Inf): incorrect ternary value %d\n", inex);
      exit (1);
    }

  /* cbrt(+/-0) =  +/-0 */
  mpfr_set_ui (x, 0, MPFR_RNDN);
  inex = mpfr_cbrt (y, x, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_NEG (y))
    {
      printf ("Error: cbrt(+0) <> +0\n");
      exit (1);
    }
  if (inex != 0)
    {
      printf ("Error: cbrt(+0): incorrect ternary value %d\n", inex);
      exit (1);
    }
  mpfr_neg (x, x, MPFR_RNDN);
  inex = mpfr_cbrt (y, x, MPFR_RNDN);
  if (MPFR_NOTZERO (y) || MPFR_IS_POS (y))
    {
      printf ("Error: cbrt(-0) <> -0\n");
      exit (1);
    }
  if (inex != 0)
    {
      printf ("Error: cbrt(-0): incorrect ternary value %d\n", inex);
      exit (1);
    }

  mpfr_set_prec (x, 53);
  mpfr_set_str (x, "8.39005285514734966412e-01", 10, MPFR_RNDN);
  mpfr_cbrt (x, x, MPFR_RNDN);
  if (mpfr_cmp_str1 (x, "9.43166207799662426048e-01"))
    {
      printf ("Error in crbrt (1)\n");
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);
  mpfr_set_str_binary (x, "0.10000100001100101001001001011001");
  mpfr_cbrt (x, x, MPFR_RNDN);
  mpfr_set_str_binary (y, "0.11001101011000100111000111111001");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in cbrt (2)\n");
      exit (1);
    }

  mpfr_set_prec (x, 32);
  mpfr_set_prec (y, 32);
  mpfr_set_str_binary (x, "-0.1100001110110000010101011001011");
  mpfr_cbrt (x, x, MPFR_RNDD);
  mpfr_set_str_binary (y, "-0.11101010000100100101000101011001");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in cbrt (3)\n");
      exit (1);
    }

  mpfr_set_prec (x, 82);
  mpfr_set_prec (y, 27);
  mpfr_set_str_binary (x, "0.1010001111011101011011000111001011001101100011110110010011011011011010011001100101e-7");
  mpfr_cbrt (y, x, MPFR_RNDD);
  mpfr_set_str_binary (x, "0.101011110001110001000100011E-2");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in cbrt (4)\n");
      exit (1);
    }

  mpfr_set_prec (x, 204);
  mpfr_set_prec (y, 38);
  mpfr_set_str_binary (x, "0.101000000001101000000001100111111011111001110110100001111000100110100111001101100111110001110001011011010110010011100101111001111100001010010100111011101100000011011000101100010000000011000101001010001001E-5");
  mpfr_cbrt (y, x, MPFR_RNDD);
  mpfr_set_str_binary (x, "0.10001001111010011011101000010110110010E-1");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in cbrt (5)\n");
      exit (1);
    }

  /* Bug (in the compiler?) found on Linux/m68k with gcc 4.0.2 */
  mpfr_set_prec (x, 5);
  mpfr_set_prec (y, 5);
  mpfr_set_str_binary (x, "1.1000E-2");
  mpfr_cbrt (y, x, MPFR_RNDN);
  mpfr_set_str_binary (x, "1.0111E-1");
  if (mpfr_cmp (x, y))
    {
      printf ("Error in cbrt (6)\n");
      exit (1);
    }

  mpfr_clear (x);
  mpfr_clear (y);
}

#define TEST_FUNCTION mpfr_cbrt
#include "tgeneric.c"

static int
cube (mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd_mode)
{
  return mpfr_pow_ui (a, b, 3, rnd_mode);
}

int
main (void)
{
  mpfr_t x;
  int r, inex;
  mpfr_prec_t p;

  tests_start_mpfr ();

  special ();

  mpfr_init (x);

  for (p = 2; p < 100; p++)
    {
      mpfr_set_prec (x, p);
      RND_LOOP (r)
        {
          mpfr_set_ui (x, 1, MPFR_RNDN);
          inex = mpfr_cbrt (x, x, (mpfr_rnd_t) r);
          if (mpfr_cmp_ui (x, 1) || inex != 0)
            {
              printf ("Error in mpfr_cbrt for x=1, rnd=%s\ngot ",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r));
              mpfr_dump (x);
              printf ("inex = %d\n", inex);
              exit (1);
            }

          mpfr_set_si (x, -1, MPFR_RNDN);
          inex = mpfr_cbrt (x, x, (mpfr_rnd_t) r);
          if (mpfr_cmp_si (x, -1) || inex != 0)
            {
              printf ("Error in mpfr_cbrt for x=-1, rnd=%s\ngot ",
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r));
              mpfr_dump (x);
              printf ("inex = %d\n", inex);
              exit (1);
            }

          if (p >= 5)
            {
              int i;
              for (i = -12; i <= 12; i++)
                {
                  mpfr_set_ui (x, 27, MPFR_RNDN);
                  mpfr_mul_2si (x, x, 3*i, MPFR_RNDN);
                  inex = mpfr_cbrt (x, x, MPFR_RNDN);
                  if (mpfr_cmp_si_2exp (x, 3, i) || inex != 0)
                    {
                      printf ("Error in mpfr_cbrt for "
                              "x = 27.0 * 2^(%d), rnd=%s\ngot ",
                              3*i, mpfr_print_rnd_mode ((mpfr_rnd_t) r));
                      mpfr_dump (x);
                      printf (", expected 3 * 2^(%d)\n", i);
                      printf ("inex = %d, expected 0\n", inex);
                      exit (1);
                    }
                }
            }
        }
    }
  mpfr_clear (x);

  test_generic (MPFR_PREC_MIN, 200, 10);

  data_check ("data/cbrt", mpfr_cbrt, "mpfr_cbrt");
  bad_cases (mpfr_cbrt, cube, "mpfr_cbrt", 256, -256, 255, 4, 128, 200, 50);

  tests_end_mpfr ();
  return 0;
}
