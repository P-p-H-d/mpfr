/* hermite -- Compute the nth degree (physicist's) Hermite polynomial.

Copyright 2025-2026 Free Software Foundation, Inc.
Contributed by Matteo Nicoli, reviewed by Paul Zimmermann.

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

/* max (x, y, z) */
#define MAX3(x,y,z) (MAX (x, MAX (y, z)))

/* extra bits used as a threshold by the small-|x| asymptotic branch */
#define MPFR_HERMITE_SMALL_X_GUARD 64

static int
asymptotic_small_x (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode,
                    mpfr_exp_t err, mpfr_prec_t realprec)
{
  mpfr_t v;
  long m, j;
  int inex_round;

  mpfr_init2 (v, realprec);

  if ((n & 1) == 0)
    {
      /* even n = 2m, m = n/2: c0 = H_n(0),
         c0(0) = 1, c0(j) = c0(j-1) * (-2(2j-1)).
         See algorithms.tex for details. */
      mpfr_set_ui (v, 1, MPFR_RNDN); /* exact */
      m = n / 2;
      for (j = 1; j <= m; j++)
         mpfr_mul_si (v, v, -2 * (2 * j - 1), MPFR_RNDN);
    }
  else
    {
      /* odd n = 2m+1, m = (n-1)/2: c1 = H'_n(0),
         c1(0) = 2, c1(j) = c1(j-1) * (-2(2j+1)),
         then lead = c1*x. See algorithms.tex for details. */
      mpfr_set_ui (v, 2, MPFR_RNDN); /* exact */
      m = (n - 1) / 2;
      for (j = 1; j <= m; j++)
        mpfr_mul_si (v, v, -2 * (2 * j + 1), MPFR_RNDN);
      mpfr_mul (v, v, x, MPFR_RNDN);
    }


  inex_round = mpfr_round_near_x (res, v, (mpfr_uexp_t) (err - 2),
                                  0, rnd_mode);

  mpfr_clear (v);

  return inex_round;
}

static int
overflow_sign (mpfr_srcptr x, long n)
{
  return (n & 1) ? MPFR_SIGN (x) : MPFR_SIGN_POS;
}

int
mpfr_hermite (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  long i;
  int x_is_zero, ternary_value = 0, inex;
  mpfr_t p1, p2, pn, first_term, second_term, f;
  mpfr_prec_t res_prec, realprec, guard_bits;
  mpfr_exp_t lost_bits;
  mpfr_exp_t b_i, f_i, g_i, h_i, q_i, a_i;

  /* these variables are used (and consequently initialized) only in the
     "Asymptotic expansion for small |x|" branch */
  mpfr_exp_t ex, l2n, rho, err;
  int is_x_tiny = 0;

  MPFR_GROUP_DECL (group);
  MPFR_SAVE_EXPO_DECL (expo);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("hermite[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  MPFR_ASSERTN(n >= 0); /* check n is non-negative */

  res_prec = MPFR_PREC (res);
  x_is_zero = MPFR_IS_ZERO (x);

  /* NaN are checke *before* any other check, according to C++ specs:
     "If the argument is NaN, NaN is returned [...]". */
  if (MPFR_IS_NAN (x))
    {
      MPFR_SET_NAN (res);
      /* as specified in the documentation, "[...] a NaN result
         (Not-a-Number) always corresponds to an exact return value." */
      MPFR_RET_NAN;
    }

  /* H_0(x) = 1. In this case, since the output is const and does not depend
     on the value of x, no further analysis on the value of x is performed */
  if (n == 0)
    {
      mpfr_set_ui (res, 1, rnd_mode);
      /* 1 is exactly representable in MPFR regardless of precision,
        so this will always return 0 */
      MPFR_RET (0);
    }

  /* we check +/-Inf after n=0 since n=0 is special */
  if (MPFR_IS_INF (x))
    {
      MPFR_SET_INF (res);
      if ((n & 1) && MPFR_IS_NEG(x))
        MPFR_SET_NEG(res);
      else
        MPFR_SET_POS(res);
      MPFR_RET (0);
    }

  /* H_n(0) when n is an odd number is always 0 */
  if (x_is_zero && (n & 1))
    {
      MPFR_SET_ZERO (res);

      /* Rationale: the coefficient of x in H_n(x) for odd n = 2m+1 is
         (-1)^m * 2 * n! / m!, which is positive for n mod 4 == 1 (m even)
         and negative for n mod 4 == 3 (m odd).
         Thus mpfr_hermite evaluates H_n(+/-0) = +/-0 for n mod 4 == 1, and
         H_n(+/-0) = -/+0 for n mod 4 == 3. */
      if (MPFR_IS_NEG(x) ^ ((n%4) == 3))
        MPFR_SET_NEG (res);
      else
        MPFR_SET_POS (res);

      MPFR_RET (0);
    }

  /* H_1(x) = 2x */
  if (n == 1)
    {
      /* result is set to 2x. The ternary value of mpfr_set is returned */
      return mpfr_mul_ui (res, x, 2, rnd_mode);
    }

  /* Taylor expansion for small |x| and n >= 2.
     For the proof of the following bound, see algorithms.tex.
     Let t = n*x^2; the tail after the leading term is bounded by the
     following geometric series:
        |tail| <= |lead| * t / (1-t),
     where |lead| is c0 (for n even), or c1*x (for n odd). */
  if (!x_is_zero && n >= 2)
    {
      /* ex = MPFR_GET_EXP(x), such that 2^(ex-1) <= |x| < 2^ex;
         l2n = ceil(log2(n)), so n <= 2^l2n;
         thus, t = n*x^2 < 2^l2n * (2^ex)^2 = 2^{l2n+2*ex}.
         We define rho = l2n+2*ex, therefore t < 2^{rho}.
         The Taylor expansion only applies for small |x|, i.e. ex < 0.
         We require ex <= -2 before computing rho: since the expansion needs
         rho <= -2 (and in fact rho very negative), any x with ex >= -1 gives
         rho >= l2n-2 >= -1 and would be rejected anyway. Requiring ex <= -2 also
         ensures that 2*ex (hence rho and err) does not overflow, since
         2*MPFR_EMIN_MIN is representable in an mpfr_exp_t whereas 2*ex for
         a large positive ex (e.g., ex close to MPFR_EMAX_MAX) would not. */
      ex = MPFR_GET_EXP (x);
      if (ex <= -2)
        {
          l2n = (mpfr_exp_t) MPFR_INT_CEIL_LOG2 (n);
          rho = l2n + 2 * ex;

          /* the bound err = -rho - 1 requires rho <= -2. In practice, require
             64 extra bits so the first rounding test usually succeeds */
          if (rho <= -2)
            {
              /* see algorithms.tex for the calculation of this error bound */
              err = -rho - 1;

              if (err >= (mpfr_exp_t) res_prec + MPFR_HERMITE_SMALL_X_GUARD)
                is_x_tiny = 1;
            }
        }
    }

  /* compute in the extended exponent range so that the final overflow/underflow
     with respect to the caller's range is handled uniformly by mpfr_check_range */
  MPFR_SAVE_EXPO_MARK (expo);

  /* analyzing all the test cases where the result is not exact (inex != 0),
     we find that the average number of bits lost per iteration, i.e.,
     lost_bits/(n-1), is about 3.27, but up to about 5.5 for n >= 20.
     We thus add 4*n guard bits for n < 20, and 6*n for n >= 20.
     For revision see eb17cda, where we have a total of 8722 such tests.
     With guard_bits * n + 10, we get a probability of failure of 0.2% */
  guard_bits = n < 20 ? 4 : 6;
  realprec = res_prec + guard_bits * n + 10;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  MPFR_GROUP_INIT_6 (group, realprec,
                     p1, p2, pn, first_term, second_term, f);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      MPFR_BLOCK_DECL (flags);

      if (is_x_tiny)
        {
          ternary_value = asymptotic_small_x (res, n, x, rnd_mode,
                                              err, realprec);
          /* if asymptotic_small_x returns 0, then it cannot round. In that
             case, our asymptotic expansion failed, so we fall back to the usual
             Ziv loop. Otherwise, we skip the main iteration and we return
             directly the result */
          if (ternary_value)
            goto clean;
        }

      i = 1;

      MPFR_BLOCK (flags, inex = mpfr_mul_ui (p1, x, 2, MPFR_RNDN));
      if (MPFR_OVERFLOW (flags))
        {
          /* the sign of H_n(x) for large |x| is that of its leading term
             (2x)^n */
          ternary_value = mpfr_overflow (res, rnd_mode,
                                         overflow_sign (x, n));
          MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
          break;
        }
      mpfr_set (f, p1, MPFR_RNDN);              /* exact */
      mpfr_set_ui (p2, 1, MPFR_RNDN);           /* exact */

      /* In the loop:
           2^a_i is a bound on the absolute error on p1 (a from algorithms.tex)
           2^b_i is a bound on the absolute error on p2 (b from algorithms.tex).
         a_i and b_i come from the previous iteration, and are initialized
         below for the first iteration (i = 1). */

      /* 2^b_i is the absolute error on p2. We use MPFR_EXP_MIN as the "minus
         infinity" exponent standing for an exact (zero) value */
      b_i = MPFR_EXP_MIN;
      /* 2^a_i is the absolute error on p1; when x = 0, p1 = 2x = 0 is exact */
      a_i = x_is_zero ? MPFR_EXP_MIN : MPFR_GET_EXP (p1) - realprec - 1;
      f_i = a_i;

      while (i < n)
        {
          /* Invariants: p1 approximates H_i(x), p2 approximates H_{i-1}(x),
             first_term will approximate 2*x*H_i(x), second_term will
             approximate 2*i*H_{i-1}(x), and pn will approximate H_{i+1}(x) */

          if (x_is_zero)
            {
              /* x = 0: the trailing term 2x*p1 is exactly 0. When i is even,
                 p2 = 0 too, so second_term = 0 and pn = 0 (all exact); when i
                 is odd, pn = -second_term = -2i*p2. In both cases we do not
                 call MPFR_GET_EXP on zero, so we use MPFR_EXP_MIN instead */
              if ((i & 1) == 0)
                {
                  MPFR_ASSERTD (MPFR_IS_ZERO (p2));
                  MPFR_SET_ZERO (pn); /* exact */
                  q_i = MPFR_EXP_MIN;
                  goto end_of_loop;
                }
              MPFR_SET_ZERO (first_term); /* first_term = 2x*p1 = 0, exact */
              h_i = MPFR_EXP_MIN;
            }
          else
            {

              /* first_term = f * p1, with absolute error at step i bounded by
                 h_i <= 2 + max(exp(first_term)-p-1, f_i+exp(p1),
                                1+MPFR_GET_EXP(x)+a_i) */
              MPFR_BLOCK (flags,
                          inex |= mpfr_mul (first_term, f, p1, MPFR_RNDN));
              if (MPFR_OVERFLOW (flags))
                break;
              h_i = 2 + MAX3 (MPFR_GET_EXP (first_term) - realprec - 1,
                              f_i + MPFR_GET_EXP (p1),
                              1 + MPFR_GET_EXP (x) + a_i);
            }

          /* second_term = p2 * 2i, with absolute error at step i bounded by
             g_i <= max(exp(second_term)-p,
                        b_i + MPFR_INT_CEIL_LOG2(2*i)+1) */
          MPFR_BLOCK (flags,
                      inex |= mpfr_mul_ui (second_term, p2, 2 * i, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          g_i = MAX (MPFR_GET_EXP (second_term) - realprec,
                     b_i + MPFR_INT_CEIL_LOG2 (2*i) + 1);

          /* pn = first_term - second_term, with absolute error at step i
             bounded by
             q_i <= 2 + max(exp(pn)-p-1, h_i, g_i)
             Note: mpfr_sub can overflow when first_term and second_term
             have opposite signs */
          MPFR_BLOCK (flags,
                      inex |= mpfr_sub (pn, first_term, second_term, MPFR_RNDN));
          if (MPFR_OVERFLOW (flags))
            break;
          q_i = 2 + MAX3 (MPFR_GET_EXP (pn) - realprec - 1, h_i, g_i);

        end_of_loop:
          /* p2 = p1, p1 = pn */
          mpfr_swap (p2, p1); /* now p2 approximates H_{i}(x) */
          mpfr_swap (p1, pn); /* now p1 approximates H_{i+1}(x) */
          b_i = a_i;          /* 2^b_i is a bound on the absolute error on p2 */
          a_i = q_i;          /* 2^a_i is a bound on the absolute error on p1 */

          i++;
        }

      /* if an overflow occurred in the recurrence (detected via flags).
         FIXME: currently, MPFR stops the iteration as soon as it encounters an
         overflow. A more refined strategy could be devised to handle intermediate
         overflows within the oscillatory region (see algorithms.tex).
         The sign is that of the leading term (2x)^n */
      if (MPFR_OVERFLOW (flags))
        {
          ternary_value = mpfr_overflow (res, rnd_mode,
                                         overflow_sign (x, n));
          MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, MPFR_FLAGS_OVERFLOW);
          break;
        }

      /* now p1 approximates H_n(x), and 2^a_i is a bound on its absolute
         error. Since ulp(p1) = 2^(EXP(p1)-realprec), we get the relative
         error is bounded by 2^(a_i - (EXP(p1) - realprec - 1)). */
      lost_bits = a_i - (MPFR_GET_EXP (p1) - realprec);

      /* if inex=0, then all the computation was exact, thus p1 is exactly
         H_n(x), otherwise we call MPFR_CAN_ROUND() to check if we can
         deduce the correct rounding */
      if (inex == 0 ||
          (lost_bits < realprec &&
            MPFR_CAN_ROUND (p1, realprec - lost_bits, res_prec, rnd_mode)))
        {
          ternary_value = mpfr_set (res, p1, rnd_mode);
          break;
        }

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_6 (group, realprec,
                           p1, p2, pn, first_term, second_term, f);
    }

 clean:
  MPFR_ZIV_FREE (loop);

  MPFR_GROUP_CLEAR (group);
  MPFR_SAVE_EXPO_FREE (expo);

  return mpfr_check_range (res, ternary_value, rnd_mode);
}
