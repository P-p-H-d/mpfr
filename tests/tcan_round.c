/* Test file for mpfr_can_round and mpfr_round_p.

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

/* Simple cases where err - prec is large enough.
   One can get failures as a consequence of r9883. */
static void
test_simple (void)
{
  int t[4] = { 2, 3, -2, -3 };  /* test powers of 2 and non-powers of 2 */
  int i;
  int r1, r2;

  for (i = 0; i < 4; i++)
    RND_LOOP (r1)
      RND_LOOP (r2)
        {
          mpfr_t b;
          int p, err, prec, inex;
          int s1, s2;
          int expected, got;

          p = 12 + (randlimb() % (2 * GMP_NUMB_BITS));
          err = p - 3;
          prec = 4;
          mpfr_init2 (b, p);
          inex = mpfr_set_si (b, t[i], MPFR_RNDN);
          MPFR_ASSERTN (inex == 0);
          got = mpfr_can_round (b, err, (mpfr_rnd_t) r1, (mpfr_rnd_t) r2, prec);
          s1 = r1;
          s2 = r2;
          if (s1 == MPFR_RNDD)
            s1 = (t[i] > 0) ? MPFR_RNDZ : MPFR_RNDA;
          if (s1 == MPFR_RNDU)
            s1 = (t[i] < 0) ? MPFR_RNDZ : MPFR_RNDA;
          if (s1 == MPFR_RNDF)
            s1 = MPFR_RNDN; /* For s1, RNDF is equivalent to RNDN. */
          if (s2 == MPFR_RNDD)
            s2 = (t[i] > 0) ? MPFR_RNDZ : MPFR_RNDA;
          if (s2 == MPFR_RNDU)
            s2 = (t[i] < 0) ? MPFR_RNDZ : MPFR_RNDA;
          /* If s1 == s2, we can round.
             s1      s2      can round
             xxx     xxx     yes
             RNDZ    RNDA    no
             RNDZ    RNDN    yes
             RNDA    RNDZ    no
             RNDA    RNDN    yes
             RNDN    RNDZ    no
             RNDN    RNDA    no
             xxx     RNDF    yes
          */
          expected = 1;
          if ((s1 == MPFR_RNDZ && s2 == MPFR_RNDA) ||
              (s1 == MPFR_RNDA && s2 == MPFR_RNDZ) ||
              (s1 == MPFR_RNDN && s2 == MPFR_RNDZ) ||
              (s1 == MPFR_RNDN && s2 == MPFR_RNDA))
            expected = 0;
          if (!!got != !!expected)
            {
              printf ("Error in test_simple for i=%d,"
                      " err=%d r1=%s, r2=%s, p=%d, prec=%d\n", i, err,
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r1),
                      mpfr_print_rnd_mode ((mpfr_rnd_t) r2), p, prec);
              printf ("b="); mpfr_dump (b);
              printf ("expected %d, got %d\n", expected, got);
              exit (1);
            }
          mpfr_clear (b);
        }
}

#define MAX_LIMB_SIZE 100

static void
check_round_p (void)
{
  mp_limb_t buf[MAX_LIMB_SIZE];
  mp_size_t n, i;
  mpfr_prec_t p;
  mpfr_exp_t err;
  int r1, r2;

  for (n = 2 ; n <= MAX_LIMB_SIZE ; n++)
    {
      /* avoid mpn_random which leaks memory */
      for (i = 0; i < n; i++)
        buf[i] = randlimb ();
      /* force the number to be normalized */
      buf[n - 1] |= MPFR_LIMB_HIGHBIT;
      p = randlimb() % ((n-1) * GMP_NUMB_BITS) + MPFR_PREC_MIN;
      err = p + randlimb () % GMP_NUMB_BITS;
      r1 = mpfr_round_p (buf, n, err, p);
      r2 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDN, MPFR_RNDZ, p);
      if (r1 != r2)
        {
          printf ("mpfr_round_p(%d) != mpfr_can_round(%d,RNDZ)!\n"
                  "bn = %ld, err0 = %ld, prec = %lu\n",
                  r1, r2, (long) n, (long) err, (unsigned long) p);
          n_trace ("b", buf, n);
          exit (1);
        }
      /* Same with RNDF: with rnd1=RNDN, rnd2=RNDF is converted to RNDN. */
      r1 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDN, MPFR_RNDN, p);
      r2 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDN, MPFR_RNDF, p);
      if (r1 != r2)
        {
          printf ("mpfr_can_round(%d,RNDN) != mpfr_can_round(%d,RNDF)!\n"
                  "bn = %ld, err0 = %ld, prec = %lu\n",
                  r1, r2, (long) n, (long) err, (unsigned long) p);
          n_trace ("b", buf, n);
          exit (1);
        }
      /* PZ: disabled those tests for now, since when {buf, n} is exactly
         representable in the target precision p, then mpfr_can_round_raw(RNDA)
         should give 0, and mpfr_can_round_raw(MPFR_RNDF) should give 1 if the
         error is small enough. */
#if 0
      /* Same with RNDF: with rnd1=RNDZ, rnd2=RNDF is converted to RNDA. */
      r1 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDZ, MPFR_RNDA, p);
      r2 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDZ, MPFR_RNDF, p);
      if (r1 != r2)
        {
          printf ("mpfr_can_round(%d,RNDA) != mpfr_can_round(%d,RNDF)!\n"
                  "bn = %ld, err0 = %ld, prec = %lu\n",
                  r1, r2, (long) n, (long) err, (unsigned long) p);
          n_trace ("b", buf, n);
          exit (1);
        }
      /* Same with RNDF: with rnd1=RNDA, rnd2=RNDF is converted to RNDZ. */
      r1 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDA, MPFR_RNDZ, p);
      r2 = mpfr_can_round_raw (buf, n, MPFR_SIGN_POS, err,
                               MPFR_RNDA, MPFR_RNDF, p);
      if (r1 != r2)
        {
          printf ("mpfr_can_round(%d,RNDZ) != mpfr_can_round(%d,RNDF)!\n"
                  "bn = %ld, err0 = %ld, prec = %lu\n",
                  r1, r2, (long) n, (long) err, (unsigned long) p);
          n_trace ("b", buf, n);
          exit (1);
        }
#endif
    }
}

/* check x=2^i with precision px, error at most 1, and target precision prec */
static void
test_pow2 (mpfr_exp_t i, mpfr_prec_t px, mpfr_rnd_t r1, mpfr_rnd_t r2,
           mpfr_prec_t prec)
{
  mpfr_t x;
  int b, expected_b, b2;

  mpfr_init2 (x, px);
  mpfr_set_ui_2exp (x, 1, i, MPFR_RNDN);
  /* for mpfr_can_round, r1=RNDF is equivalent to r1=RNDN (the sign of the
     error is not known) */
  b = !!mpfr_can_round (x, i+1, r1, r2, prec);
  /* Note: If mpfr_can_round succeeds for both
     (r1,r2) = (MPFR_RNDD,MPFR_RNDN) and
     (r1,r2) = (MPFR_RNDU,MPFR_RNDN), then it should succeed for
     (r1,r2) = (MPFR_RNDN,MPFR_RNDN). So, the condition on prec below
     for r1 = MPFR_RNDN should be the most restrictive between those
     for r1 = any directed rounding mode.
     For r1 like MPFR_RNDA, the unrounded, unknown number may be anyone
     in [2^i-1,i]. As both 2^i-1 and 2^i fit on i bits, one cannot round
     in any precision >= i bits, hence the condition prec < i; prec = i-1
     will work here for r2 = MPFR_RNDN thanks to the even-rounding rule
     (and also with rounding ties away from zero). */
  expected_b =
    MPFR_IS_LIKE_RNDD (r1, MPFR_SIGN_POS) ?
    (MPFR_IS_LIKE_RNDU (r2, MPFR_SIGN_POS) ? 0 : prec <= i) :
    MPFR_IS_LIKE_RNDU (r1, MPFR_SIGN_POS) ?
    (MPFR_IS_LIKE_RNDD (r2, MPFR_SIGN_POS) ? 0 : prec < i) :
    (r2 != MPFR_RNDN ? 0 : prec < i);
  if (b != expected_b)
    {
      printf ("Error for x=2^%d, px=%lu, err=%d, r1=%s, r2=%s, prec=%d\n",
              (int) i, (unsigned long) px, (int) i + 1,
              mpfr_print_rnd_mode (r1), mpfr_print_rnd_mode (r2), (int) prec);
      printf ("Expected %d, got %d\n", expected_b, b);
      exit (1);
    }

  if (r1 == MPFR_RNDN && r2 == MPFR_RNDZ)
    {
      /* Similar test to the one done in src/round_p.c
         for MPFR_WANT_ASSERT >= 2. */
      b2 = !!mpfr_round_p (MPFR_MANT(x), MPFR_LIMB_SIZE(x), i+1, prec);
      if (b2 != b)
        {
          printf ("Error for x=2^%d, px=%lu, err=%d, prec=%d\n",
                  (int) i, (unsigned long) px, (int) i + 1, (int) prec);
          printf ("mpfr_can_round gave %d, mpfr_round_p gave %d\n", b, b2);
          exit (1);
        }
    }

  mpfr_clear (x);
}

static void
check_can_round (void)
{
  mpfr_t x, xinf, xsup, yinf, ysup;
  int precx, precy, err;
  int rnd1, rnd2;
  int i, u[3] = { 0, 1, 256 };
  int inex;
  int expected, got;
  int maxerr;

  mpfr_inits2 (4 * GMP_NUMB_BITS, x, xinf, xsup, yinf, ysup, (mpfr_ptr) 0);

  for (precx = 3 * GMP_NUMB_BITS - 3; precx <= 3 * GMP_NUMB_BITS + 3; precx++)
    {
      mpfr_set_prec (x, precx);
      for (precy = precx - 4; precy <= precx + 4; precy++)
        {
          mpfr_set_prec (yinf, precy);
          mpfr_set_prec (ysup, precy);

          for (i = 0; i <= 3; i++)
            {
              if (i <= 2)
                {
                  /* Test x = 2^(precx - 1) + u */
                  mpfr_set_ui_2exp (x, 1, precx - 1, MPFR_RNDN);
                  mpfr_add_ui (x, x, u[i], MPFR_RNDN);
                }
              else
                {
                  /* Test x = 2^precx - 1 */
                  mpfr_set_ui_2exp (x, 1, precx, MPFR_RNDN);
                  mpfr_sub_ui (x, x, 1, MPFR_RNDN);
                }
              MPFR_ASSERTN (mpfr_get_exp (x) == precx);

              maxerr = precy + 3;
              if (4 * GMP_NUMB_BITS < maxerr)
                maxerr = 4 * GMP_NUMB_BITS;
              for (err = precy; err <= maxerr; err++)
                {
                  mpfr_set_ui_2exp (xinf, 1, precx - err, MPFR_RNDN);
                  inex = mpfr_add (xsup, x, xinf, MPFR_RNDN);
                  /* Since EXP(x) = precx, and xinf = 2^(precx-err),
                     x + xinf is exactly representable on 4 * GMP_NUMB_BITS
                     nbits as long as err <= 4 * GMP_NUMB_BITS */
                  MPFR_ASSERTN (inex == 0);
                  inex = mpfr_sub (xinf, x, xinf, MPFR_RNDN);
                  MPFR_ASSERTN (inex == 0);
                  RND_LOOP (rnd1)
                    RND_LOOP (rnd2)
                      {
                        /* TODO: Test r2 == MPFR_RNDF. The following "continue"
                           was added while this case had not been specified
                           yet, but this is no longer the case. */
                        if (rnd2 == MPFR_RNDF)
                          continue;
                        mpfr_set (yinf, MPFR_IS_LIKE_RNDD (rnd1, 1) ?
                                  x : xinf, (mpfr_rnd_t) rnd2);
                        mpfr_set (ysup, MPFR_IS_LIKE_RNDU (rnd1, 1) ?
                                  x : xsup, (mpfr_rnd_t) rnd2);
                        expected = !! mpfr_equal_p (yinf, ysup);
                        got = !! mpfr_can_round (x, err, (mpfr_rnd_t) rnd1,
                                                 (mpfr_rnd_t) rnd2, precy);
                        if (got != expected)
                          {
                            printf ("Error in check_can_round on:\n"
                                    "precx=%d, precy=%d, i=%d, err=%d, "
                                    "rnd1=%s, rnd2=%s: expected %d, got %d\n",
                                    precx, precy, i, err,
                                    mpfr_print_rnd_mode ((mpfr_rnd_t) rnd1),
                                    mpfr_print_rnd_mode ((mpfr_rnd_t) rnd2),
                                    expected, got);
                            printf ("x="); mpfr_dump (x);
                            printf ("yinf="); mpfr_dump (yinf);
                            printf ("ysup="); mpfr_dump (ysup);
                            exit (1);
                          }
                      }
                }
            }
        }
    }

  mpfr_clears (x, xinf, xsup, yinf, ysup, (mpfr_ptr) 0);
}

/* test of RNDNA (nearest with ties to away) */
static void
test_rndna (void)
{
  mpfr_t x;
  int inex;

  mpfr_init2 (x, 10);
  mpfr_set_str_binary (x, "1111111101"); /* 1021 */
  inex = mpfr_prec_round (x, 9, MPFR_RNDNA);
  MPFR_ASSERTN(inex > 0);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 1022) == 0);
  mpfr_set_prec (x, 10);
  mpfr_set_str_binary (x, "1111111101"); /* 1021 */
  inex = mpfr_prec_round (x, 9, MPFR_RNDN);
  MPFR_ASSERTN(inex < 0);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 1020) == 0);
  mpfr_set_prec (x, 10);
  mpfr_set_str_binary (x, "1111111011"); /* 1019 */
  inex = mpfr_prec_round (x, 9, MPFR_RNDNA);
  MPFR_ASSERTN(inex > 0);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 1020) == 0);
  mpfr_set_prec (x, 10);
  mpfr_set_str_binary (x, "1111111011"); /* 1019 */
  inex = mpfr_prec_round (x, 9, MPFR_RNDN);
  MPFR_ASSERTN(inex > 0);
  MPFR_ASSERTN(mpfr_cmp_ui (x, 1020) == 0);
  mpfr_clear (x);
}

int
main (void)
{
  mpfr_t x;
  mpfr_prec_t i, j, k;
  int r1, r2;
  int n;

  tests_start_mpfr ();

  test_rndna ();
  test_simple ();

  /* checks that rounds to nearest sets the last
     bit to zero in case of equal distance */
  mpfr_init2 (x, 59);
  mpfr_set_str_binary (x, "-0.10010001010111000011110010111010111110000000111101100111111E663");
  if (mpfr_can_round (x, 54, MPFR_RNDZ, MPFR_RNDZ, 53))
    {
      printf ("Error (1) in mpfr_can_round\n");
      exit (1);
    }

  mpfr_set_str_binary (x, "-Inf");
  if (mpfr_can_round (x, 2000, MPFR_RNDZ, MPFR_RNDZ, 2000))
    {
      printf ("Error (2) in mpfr_can_round\n");
      exit (1);
    }

  mpfr_set_prec (x, 64);
  mpfr_set_str_binary (x, "0.1011001000011110000110000110001111101011000010001110011000000000");
  if (mpfr_can_round (x, 65, MPFR_RNDN, MPFR_RNDN, 54))
    {
      printf ("Error (3) in mpfr_can_round\n");
      exit (1);
    }

  mpfr_set_prec (x, 137);
  mpfr_set_str_binary (x, "-0.10111001101001010110011000110100111010011101101010010100101100001110000100111111011101010110001010111100100101110111100001000010000000000E-97");
  if (! mpfr_can_round (x, 132, MPFR_RNDU, MPFR_RNDZ, 128))
    {
      printf ("Error (4) in mpfr_can_round\n");
      exit (1);
    }

  /* in the following, we can round but cannot determine the inexact flag */
  mpfr_set_prec (x, 86);
  mpfr_set_str_binary (x, "-0.11100100010011001111011010100111101010011000000000000000000000000000000000000000000000E-80");
  if (! mpfr_can_round (x, 81, MPFR_RNDU, MPFR_RNDZ, 44))
    {
      printf ("Error (5) in mpfr_can_round\n");
      exit (1);
    }

  mpfr_set_prec (x, 62);
  mpfr_set_str (x, "0.ff4ca619c76ba69", 16, MPFR_RNDZ);
  for (i = 30; i < 99; i++)
    for (j = 30; j < 99; j++)
      RND_LOOP (r1)
        RND_LOOP (r2)
          {
            /* test for assertions */
            mpfr_can_round (x, i, (mpfr_rnd_t) r1, (mpfr_rnd_t) r2, j);
          }

  test_pow2 (32, 32, MPFR_RNDN, MPFR_RNDN, 32);
  test_pow2 (174, 174, MPFR_RNDN, MPFR_RNDN, 174);
  test_pow2 (174, 174, MPFR_RNDU, MPFR_RNDN, 174);
  test_pow2 (176, 129, MPFR_RNDU, MPFR_RNDU, 174);
  test_pow2 (176, 2, MPFR_RNDZ, MPFR_RNDZ, 174);
  test_pow2 (176, 2, MPFR_RNDU, MPFR_RNDU, 176);

  /* Tests for x = 2^i (E(x) = i+1) with error at most 1 = 2^0. */
  for (n = 0; n < 100; n++)
    {
      /* TODO: Test r2 == MPFR_RNDF (add its support in test_pow2). The
         exclusion below was added while this case had not been specified
         yet, but this is no longer the case. */
      i = (randlimb() % 200) + 4;
      for (j = i - 2; j < i + 2; j++)
        RND_LOOP (r1)
          RND_LOOP (r2)
            if (r2 != MPFR_RNDF)
              for (k = MPFR_PREC_MIN; k <= i + 2; k++)
                test_pow2 (i, k, (mpfr_rnd_t) r1, (mpfr_rnd_t) r2, j);
    }

  mpfr_clear (x);

  check_can_round ();

  check_round_p ();

  tests_end_mpfr ();
  return 0;
}
