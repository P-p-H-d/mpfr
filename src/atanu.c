/* mpfr_atanu  -- atanu(x)  = atan(x)*u/(2*pi)
   mpfr_atanpi -- atanpi(x) = atan(x)/pi

Copyright 2021-2025 Free Software Foundation, Inc.
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

#define MPFR_NEED_LONGLONG_H
#include "mpfr-impl.h"

/* put in y the correctly rounded value of atan(x)*u/(2*pi) */
int
mpfr_atanu (mpfr_ptr y, mpfr_srcptr x, unsigned long u, mpfr_rnd_t rnd_mode)
{
  mpfr_t tmp, pi;
  mpfr_prec_t prec;
  mpfr_exp_t expx;
  int inex;
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg u=%lu rnd=%d", mpfr_get_prec(x), mpfr_log_prec, x, u,
      rnd_mode),
     ("y[%Pd]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y,
      inex));

  /* Singular cases */
  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (x)))
    {
      if (MPFR_IS_NAN (x))
        {
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x))
        {
          /* atanu(+Inf,u) = u/4, atanu(-Inf,u) = -u/4 */
          if (MPFR_IS_POS (x))
            return mpfr_set_ui_2exp (y, u, -2, rnd_mode);
          else
            {
              inex = mpfr_set_ui_2exp (y, u, -2, MPFR_INVERT_RND (rnd_mode));
              MPFR_CHANGE_SIGN (y);
              return -inex;
            }
        }
      else /* necessarily x=0 */
        {
          MPFR_ASSERTD(MPFR_IS_ZERO(x));
          /* atan(0)=0 with same sign, even when u=0 to ensure
             atanu(-x,u) = -atanu(x,u) */
          MPFR_SET_ZERO (y);
          MPFR_SET_SAME_SIGN (y, x);
          MPFR_RET (0); /* exact result */
        }
    }

  if (u == 0) /* return 0 with sign of x, which is coherent with case x=0 */
    {
      MPFR_SET_ZERO (y);
      MPFR_SET_SAME_SIGN (y, x);
      MPFR_RET (0);
    }

  if (mpfr_cmpabs_ui (x, 1) == 0)
    {
      /* |x| = 1: atanu(1,u) = u/8, atanu(-1,u)=-u/8 */
      /* we can't use mpfr_set_si_2exp with -u since -u might not be
         representable as long */
      if (MPFR_SIGN(x) > 0)
        return mpfr_set_ui_2exp (y, u, -3, rnd_mode);
      else
        {
          inex = mpfr_set_ui_2exp (y, u, -3, MPFR_INVERT_RND(rnd_mode));
          MPFR_CHANGE_SIGN(y);
          return -inex;
        }
    }

  /* For x>=1, we have pi/2-1/x < atan(x) < pi/2, thus
     u/4-u/(2*pi*x) < atanu(x,u) < u/4, and the relative difference between
     atanu(x,u) and u/4 is less than 2/(pi*x) < 1/x <= 2^(1-EXP(x)).
     If the relative difference is <= 2^(-prec-2), then the difference
     between atanu(x,u) and u/4 is <= 1/4*ulp(u/4) <= 1/2*ulp(RN(u/4)).
     We also require x >= 2^64, which implies x > 2*u/pi, so that
     (u-1)/4 < u/4-u/(2*pi*x) < u/4. */
  expx = MPFR_GET_EXP(x);
  if (expx >= 65 && expx - 1 >= MPFR_PREC(y) + 2)
    {
      prec = (MPFR_PREC(y) <= 63) ? 65 : MPFR_PREC(y) + 2;
      /* now prec > 64 and prec > MPFR_PREC(y)+1 */
      mpfr_init2 (tmp, prec);
      /* since expx >= 65, we have emax >= 65, thus u is representable here,
         and we don't need to work in an extended exponent range */
      inex = mpfr_set_ui (tmp, u, MPFR_RNDN); /* exact since prec >= 64 */
      MPFR_ASSERTD(inex == 0);
      mpfr_nextbelow (tmp);
      /* Since prec >= 65, the last significant bit of tmp is 1, and since
         prec > PREC(y), tmp is not representable in the target precision,
         which ensures we will get a correct ternary value below. */
      MPFR_ASSERTD(mpfr_min_prec(tmp) > MPFR_PREC(y));
      if (MPFR_SIGN(x) < 0)
        MPFR_CHANGE_SIGN(tmp);
      /* since prec >= PREC(y)+2, the rounding of tmp is correct */
      inex = mpfr_div_2ui (y, tmp, 2, rnd_mode);
      mpfr_clear (tmp);
      return inex;
    }

  prec = MPFR_PREC (y);

  MPFR_SAVE_EXPO_MARK (expo);

  prec += MPFR_INT_CEIL_LOG2(prec) + 10;

  mpfr_init2 (tmp, prec);
  mpfr_init2 (pi, prec);

  MPFR_ZIV_INIT (loop, prec);
  for (;;)
    {
      /* In the error analysis below, each thetax denotes a variable such that
         |thetax| <= 2^(1-prec) */
      mpfr_atan (tmp, x, MPFR_RNDA);
      /* tmp = atan(x) * (1 + theta1), and tmp cannot be zero since we rounded
         away from zero, and the case x=0 was treated before */
      /* first multiply by u to avoid underflow issues */
      mpfr_mul_ui (tmp, tmp, u, MPFR_RNDA);
      /* tmp = atan(x)*u * (1 + theta2)^2, and |tmp| >= 0.5*2^emin */
      mpfr_const_pi (pi, MPFR_RNDZ); /* round toward zero since we we will
                                        divide by pi, to round tmp away */
      /* pi = Pi * (1 + theta3) */
      mpfr_div (tmp, tmp, pi, MPFR_RNDA);
      /* tmp = atan(x)*u/Pi * (1 + theta4)^4, with |tmp| > 0 */
      /* since we rounded away from 0, if we get 0.5*2^emin here, it means
         |atanu(x,u)| < 0.25*2^emin (pi is not exact) thus we have underflow */
      if (MPFR_EXP(tmp) == __gmpfr_emin)
        {
          /* mpfr_underflow rounds away for RNDN */
          mpfr_clear (tmp);
          mpfr_clear (pi);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_underflow (y,
                            (rnd_mode == MPFR_RNDN) ? MPFR_RNDZ : rnd_mode, 1);
        }
      mpfr_div_2ui (tmp, tmp, 1, MPFR_RNDA); /* exact */
      /* tmp = atan(x)*u/(2*Pi) * (1 + theta4)^4 */
      /* since |(1 + theta4)^4 - 1| <= 8*|theta4| for prec >= 3,
         the relative error is less than 2^(4-prec) */
      MPFR_ASSERTD(!MPFR_IS_ZERO(tmp));
      if (MPFR_LIKELY (MPFR_CAN_ROUND (tmp, prec - 4,
                                       MPFR_PREC (y), rnd_mode)))
        break;
      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (tmp, prec);
      mpfr_set_prec (pi, prec);
    }
  MPFR_ZIV_FREE (loop);

  inex = mpfr_set (y, tmp, rnd_mode);
  mpfr_clear (tmp);
  mpfr_clear (pi);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inex, rnd_mode);
}

int
mpfr_atanpi (mpfr_ptr y, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  return mpfr_atanu (y, x, 2, rnd_mode);
}
