/* mpfr_compound_si --- compound(x,n) = (1+x)^n

Copyright 2021-2023 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

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
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#define MPFR_NEED_LONGLONG_H /* needed for MPFR_INT_CEIL_LOG2 */
#include "mpfr-impl.h"

/* assuming |(1+x)^n - 1| < 1/4*ulp(1), return correct rounding,
   where s is the sign of n*log2(1+x) */
static int
mpfr_compound_near_one (mpfr_ptr y, int s, mpfr_rnd_t rnd_mode)
{
  mpfr_set_ui (y, 1, rnd_mode); /* exact */
  if (rnd_mode == MPFR_RNDN || rnd_mode == MPFR_RNDF
      || (s > 0 && (rnd_mode == MPFR_RNDZ || rnd_mode == MPFR_RNDD))
      || (s < 0 && (rnd_mode == MPFR_RNDA || rnd_mode == MPFR_RNDU)))
    {
      /* round toward 1 */
      return -s;
    }
  else if (s > 0) /* necessarily RNDA or RNDU */
    {
      /* round toward +Inf */
      mpfr_nextabove (y);
      return +1;
    }
  else /* necessarily s < 0 and RNDZ or RNDD */
    {
      /* round toward 0 */
      mpfr_nextbelow (y);
      return -1;
    }
}

/* put in y the correctly rounded value of (1+x)^n */
int
mpfr_compound_si (mpfr_ptr y, mpfr_srcptr x, long n, mpfr_rnd_t rnd_mode)
{
  int inexact, compared, k, nloop;
  mpfr_t t, u;
  mpfr_prec_t prec, extra;
  mpfr_rnd_t rnd1;
  MPFR_ZIV_DECL (loop);
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_LOG_FUNC
    (("x[%Pu]=%.*Rg n=%ld rnd=%d",
      mpfr_get_prec(x), mpfr_log_prec, x, n, rnd_mode),
     ("y[%Pu]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y, inexact));

  /* Special cases */
  if (MPFR_IS_SINGULAR (x))
    {
      if (MPFR_IS_INF (x) && MPFR_IS_NEG (x))
        {
          /* compound(-Inf,n) is NaN */
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (n == 0 || MPFR_IS_ZERO (x))
        {
          /* compound(x,0) = 1 for x >= -1 or NaN (the only special value
             of x that is not concerned is -Inf, already handled);
             compound(0,n) = 1 */
          return mpfr_set_ui (y, 1, rnd_mode);
        }
      else if (MPFR_IS_NAN (x))
        {
          /* compound(NaN,n) is NaN, except for n = 0, already handled. */
          MPFR_SET_NAN (y);
          MPFR_RET_NAN;
        }
      else if (MPFR_IS_INF (x)) /* x = +Inf */
        {
          MPFR_ASSERTD (MPFR_IS_POS (x));
          if (n < 0) /* (1+Inf)^n = +0 for n < 0 */
            MPFR_SET_ZERO (y);
          else /* n > 0: (1+Inf)^n = +Inf */
            MPFR_SET_INF (y);
          MPFR_SET_POS (y);
          MPFR_RET (0); /* exact 0 or infinity */
        }
    }

  /* (1+x)^n = NaN for x < -1 */
  compared = mpfr_cmp_si (x, -1);
  if (compared < 0)
    {
      MPFR_SET_NAN (y);
      MPFR_RET_NAN;
    }

  /* compound(x,0) gives 1 for x >= 1 */
  if (n == 0)
    return mpfr_set_ui (y, 1, rnd_mode);

  if (compared == 0)
    {
      if (n < 0)
        {
          /* compound(-1,n) = +Inf with divide-by-zero exception */
          MPFR_SET_INF (y);
          MPFR_SET_POS (y);
          MPFR_SET_DIVBY0 ();
          MPFR_RET (0);
        }
      else
        {
          /* compound(-1,n) = +0 */
          MPFR_SET_ZERO (y);
          MPFR_SET_POS (y);
          MPFR_RET (0);
        }
    }

  if (n == 1)
    return mpfr_add_ui (y, x, 1, rnd_mode);

  MPFR_SAVE_EXPO_MARK (expo);

  prec = MPFR_PREC(y);
  prec += MPFR_INT_CEIL_LOG2 (prec) + 6;

  mpfr_init2 (t, prec);
  mpfr_init2 (u, prec);

  k = MPFR_INT_CEIL_LOG2(SAFE_ABS (unsigned long, n));  /* thus |n| <= 2^k */

  /* we compute u=log2p1(x) with prec+extra bits, since we loose some bits
     in 2^u */
  extra = 0;
  rnd1 = VSIGN (n) == MPFR_SIGN (x) ? MPFR_RNDD : MPFR_RNDU;

  MPFR_ZIV_INIT (loop, prec);
  for (nloop = 0; ; nloop++)
    {
      unsigned int inex;
      mpfr_exp_t e, e2;
      mpfr_prec_t precu = MPFR_ADD_PREC (prec, extra);
      mpfr_prec_t new_extra;
      mpfr_rnd_t rnd2;

      /* We compute (1+x)^n as 2^(n*log2p1(x)),
         and we round toward 1, thus we round n*log2p1(x) toward 0,
         thus for x*n > 0 we round log2p1(x) toward -Inf, and for x*n < 0
         we round log2p1(x) toward +Inf. */
      inex = mpfr_log2p1 (u, x, rnd1) != 0;
      e = MPFR_GET_EXP (u);
      /* |u - log2(1+x)| <= ulp(t) = 2^(e-precu) */
      inex |= mpfr_mul_si (u, u, n, MPFR_RNDZ) != 0;
      e2 = MPFR_GET_EXP (u);
      /* |u - n*log2(1+x)| <= 2^(e2-precu) + |n|*2^(e-precu)
                           <= 2^(e2-precu) + 2^(e+k-precu) <= 2^(e+k+1-precu)
                          where |n| <= 2^k, and e2 is the new exponent of u. */
      MPFR_ASSERTD (e2 <= e + k);
      e += k + 1;
      MPFR_ASSERTN (e2 <= MPFR_PREC_MAX);
      new_extra = e2 > 0 ? e2 : 0;
      /* |u - n*log2(1+x)| <= 2^(e-precu) */
      /* detect overflow: since we rounded n*log2p1(x) toward 0,
         if n*log2p1(x) >= __gmpfr_emax, we are sure there is overflow. */
      if (mpfr_cmp_si (u, __gmpfr_emax) >= 0)
        {
          MPFR_ZIV_FREE (loop);
          mpfr_clear (t);
          mpfr_clear (u);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_overflow (y, rnd_mode, 1);
        }
      /* detect underflow: similarly, since we rounded n*log2p1(x) toward 0,
         if n*log2p1(x) < __gmpfr_emin-1, we are sure there is underflow. */
      if (mpfr_cmp_si (u, __gmpfr_emin - 1) < 0)
        {
          MPFR_ZIV_FREE (loop);
          mpfr_clear (t);
          mpfr_clear (u);
          MPFR_SAVE_EXPO_FREE (expo);
          return mpfr_underflow (y,
                            rnd_mode == MPFR_RNDN ? MPFR_RNDZ : rnd_mode, 1);
        }
      /* Detect cases where result is 1 or 1+ulp(1) or 1-1/2*ulp(1):
         |2^u - 1| = |exp(u*log(2)) - 1| <= |u|*log(2) < |u| */
      if (nloop == 0 && MPFR_GET_EXP(u) < - (mpfr_exp_t) MPFR_PREC(y))
        {
          /* since ulp(1) = 2^(1-PREC(y)), we have |u| < 1/4*ulp(1) */
          /* mpfr_compound_near_one must be called in the extended
             exponent range, so that 1 is representable. */
          inexact = mpfr_compound_near_one (y, MPFR_SIGN (u), rnd_mode);
          goto end;
        }
      /* FIXME: mpfr_exp2 could underflow to the smallest positive number
         since MPFR_RNDA is used, and this case will not be detected by
         MPFR_CAN_ROUND (see BUGS). Either fix that, or do early underflow
         detection (which may be necessary). */
      /* round 2^u toward 1 */
      rnd2 = MPFR_IS_POS (u) ? MPFR_RNDD : MPFR_RNDU;
      inex |= mpfr_exp2 (t, u, rnd2) != 0;
      /* we had |u - n*log2(1+x)| < 2^(e-precu)
         thus u = n*log2(1+x) + delta with |delta| < 2^(e-precu)
         then 2^u = (1+x)^n * 2^delta with |delta| < 2^(e-precu).
         For |delta| < 0.5, |2^delta - 1| <= |delta| thus
         |t - (1+x)^n| <= ulp(t) + |t|*2^(e-precu)
                       < 2^(EXP(t)-prec) + 2^(EXP(t)+e-precu) */
      e = (precu - prec >= e) ? 1 : e + 1 - (precu - prec);
      /* now |t - (1+x)^n| < 2^(EXP(t)+e-prec) */

      if (MPFR_LIKELY (!inex ||
                       MPFR_CAN_ROUND (t, prec - e, MPFR_PREC(y), rnd_mode)))
        break;

      /* If t fits in the target precision (or with 1 more bit), then we can
         round, assuming the working precision is large enough, but the above
         MPFR_CAN_ROUND() will fail because we cannot determine the ternary
         value. However since we rounded t toward 1, we can determine it.
         Since the error in the approximation t is at most 2^e ulp(t),
         this error should be less than 1/2 ulp(y), thus we should have
         prec - PREC(y) >= e + 1. */
      if (mpfr_min_prec (t) <= MPFR_PREC(y) + 1 &&
          prec - MPFR_PREC(y) >= e + 1)
        {
          /* we add/subtract one ulp to get the correct rounding */
          if (rnd2 == MPFR_RNDD) /* t was rounded downwards */
            mpfr_nextabove (t);
          else
            mpfr_nextbelow (t);
          break;
        }

      /* Check if x^n fits in the target precision (or 1 more bit for
         rounding to nearest),
         then if x is very large Ziv's strategy will take too long.
         If x fits into k bits, then 2^(k-1) <= x < 2^k, thus
         2^(n*(k-1)) <= x^n < 2^(k*n), thus x^n has between n*(k-1)+1 and
         k*n bits. Since this does not depend on the working precision,
         we only check this at the first iteration. */
      if (nloop == 0 && n > 1 && mpfr_cmp_ui (x, 65535) > 0)
        {
          mpfr_prec_t kx = mpfr_min_prec (x);
          mpfr_prec_t p = MPFR_PREC(y) + (rnd_mode == MPFR_RNDN);
          if (n * (kx - 1) + 1 <= p)
            {
              /* first check that x^n really fits into p bits */
              mpfr_t v;
              mpfr_prec_t min_prec_v;

              mpfr_init2 (v, p);
              inexact = mpfr_pow_ui (v, x, n, MPFR_RNDZ);
              min_prec_v = mpfr_min_prec (v);
              mpfr_clear (v);
              if (inexact == 0)
                {
                  /* (x+1)^n = x^n * (1 + 1/x)^n
                     For directed rounding, we can round when (1 + 1/x)^n
                     < 1 + 2^-p, and then the result is x^n,
                     except for rounding up. Indeed, if (1 + 1/x)^n < 1 + 2^-p,
                     1 <= (x+1)^n < x^n * (1 + 2^-p) = x^n + x^n/2^p
                     < x^n + ulp(x^n).
                     For rounding to nearest, we can round when (1 + 1/x)^n
                     < 1 + 2^-p, and then the result is x^n when x^n fits
                     into p-1 bits, and nextabove(x^n) otherwise. */
                  mpfr_ui_div (t, 1, x, MPFR_RNDU);
                  mpfr_add_ui (t, t, 1, MPFR_RNDU);
                  mpfr_pow_ui (t, t, n, MPFR_RNDU);
                  mpfr_sub_ui (t, t, 1, MPFR_RNDU);
                  /* t cannot be zero */
                  if (MPFR_GET_EXP(t) < -MPFR_PREC(y))
                    {
                      mpfr_pow_si (y, x, n, MPFR_RNDZ);
                      if (rnd_mode == MPFR_RNDN && min_prec_v == p)
                        rnd_mode = MPFR_RNDU; /* midpoint: round up */
                      if (rnd_mode != MPFR_RNDU && rnd_mode != MPFR_RNDA)
                        inexact = -1;
                      else /* round up */
                        {
                          mpfr_nextabove (y);
                          inexact = 1;
                        }
                      goto end;
                    }
                }
            }
        }

      /* Exact cases like compound(0.5,2) = 9/4 must be detected, since
         except for 1+x power of 2, the log2p1 above will be inexact,
         so that in the Ziv test, inexact != 0 and MPFR_CAN_ROUND will
         fail (even for RNDN, as the ternary value cannot be determined),
         yielding an infinite loop.
         For an exact case in precision prec(y), 1+x will necessarily
         be exact in precision prec(y), thus also in prec(t), where
         prec(t) >= prec(y), and we can use mpfr_pow_si under this
         condition (which will also evaluate some non-exact cases). */
      if (mpfr_add_ui (t, x, 1, MPFR_RNDZ) == 0)
        {
          inexact = mpfr_pow_si (y, t, n, rnd_mode);
          goto end;
        }

      MPFR_ZIV_NEXT (loop, prec);
      mpfr_set_prec (t, prec);
      extra = new_extra;
      mpfr_set_prec (u, MPFR_ADD_PREC (prec, extra));
    }

  inexact = mpfr_set (y, t, rnd_mode);

 end:
  MPFR_ZIV_FREE (loop);
  mpfr_clear (t);
  mpfr_clear (u);

  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
