/* legendre -- Compute the nth degree Legendre polynomial.

Copyright 2025-2026 Free Software Foundation, Inc.
Contributed by Matteo Nicoli and Paul Zimmermann.

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
#define MPFR_LEGENDRE_SMALL_X_GUARD 64

static int
asymptotic_small_x (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode,
                    int *inex_round, mpfr_exp_t err)
{
  mpfr_t v;
  long m, j;
  unsigned inex;
  mpfr_prec_t res_prec, realprec;

  MPFR_GROUP_DECL (small_x);
  MPFR_ZIV_DECL (loop);

  res_prec = MPFR_PREC (res);

  /* a-priori rounding error bound. See algorithms.tex for details.
     After n operations the error is at most n*2*ulp(v) <= 2^(ell+1)*ulp(v),
     where ell = ceil(log2(n)); 15 extra bits has been added for safety */
  realprec = res_prec + MPFR_INT_CEIL_LOG2 (n) + 15;

  MPFR_GROUP_INIT_1 (small_x, realprec, v);
  MPFR_ZIV_INIT (loop, realprec);

  for (;;)
    {
      /* v = 1 is exact */
      inex = mpfr_set_ui (v, 1, MPFR_RNDN);

      if ((n & 1) == 0)
        {
          /* Even n = 2m, m = n/2: c0 = P_n(0),
             c0(0) = 1, c0(j) = c0(j-1) * (-(2j-1)) / (2j).
             See algorithms.tex for details */
          m = n / 2;
          for (j = 1; j <= m; j++)
            {
              inex |= mpfr_mul_si (v, v, -(2 * j - 1), MPFR_RNDN);
              inex |= mpfr_div_ui (v, v, 2 * j, MPFR_RNDN);
            }
        }
      else
        {
          /* Odd n = 2m+1, m = (n-1)/2: c1 = P'_n(0),
             c1(0) = 1, c1(j) = c1(j-1) * (-(2j+1)) / (2j),
             then lead = c1*x. See algorithms.tex for details */
          m = (n - 1) / 2;
          for (j = 1; j <= m; j++)
            {
              inex |= mpfr_mul_si (v, v, -(2 * j + 1), MPFR_RNDN);
              inex |= mpfr_div_ui (v, v, 2 * j, MPFR_RNDN);
            }
          inex |= mpfr_mul (v, v, x, MPFR_RNDN);
        }

      /* if inex=0, then all the computation was exact, thus v is exactly V,
         otherwise we call MPFR_CAN_ROUND() to check if we can deduce
         the correct rounding */
      if (!inex || MPFR_CAN_ROUND (v, realprec - err, res_prec, rnd_mode))
        break;

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_1 (small_x, realprec, v);
    }

  MPFR_ZIV_FREE (loop);

  *inex_round = mpfr_round_near_x (res, v, (mpfr_uexp_t) (err - 2),
                                   0, rnd_mode);

  MPFR_GROUP_CLEAR (small_x);

  return inex;
}

int
mpfr_legendre (mpfr_ptr res, long n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  int ternary_value = 0;
  unsigned is_within_domain = 1, inex;
  mpfr_prec_t res_prec, realprec;

  /* these variables are used (and consequently initialized) only in the
     "Asymptotic expansion for small |x|" branch */
  mpfr_exp_t ex, l2n, rho, err;
    int inex_round;

  /* the following variables are used (and consequently initialized) only
     for n >= 2, where x is not equal to -1, 0 or 1 */
  long i;
  mpfr_t p1, p2, pn, first_term, second_term;
  mpfr_exp_t lost_bits;
  mpfr_exp_t b_i, log2_i_m1, g_i, h_i, q_i, a_i, a_n;
  int x_is_zero;

  MPFR_GROUP_DECL (main);
  MPFR_ZIV_DECL (loop);

  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", MPFR_PREC (x), mpfr_log_prec, x, rnd_mode),
     ("legendre[%Pd]=%.*Rg", MPFR_PREC (res), mpfr_log_prec, res));

  MPFR_ASSERTN(n >= 0); /* check n is non-negative */

  /* first, check if x belongs to the domain [-1,1].
     If it's not, res is set to NAN, and 0 is returned */
  is_within_domain &= mpfr_lessequal_p (x, __gmpfr_one);
  is_within_domain &= mpfr_greaterequal_p (x, __gmpfr_mone);

  if (!is_within_domain)
    {
      MPFR_SET_NAN (res);
      MPFR_RET_NAN;
    }

  /* 1 and -1 are the (respectively) upper and lower bound of the Legendre
     polynomial's canonical domain. We can evaluate Pn (for any n) without
     using Bonnet's recursion */
  if (mpfr_equal_p (x, __gmpfr_one))
    return mpfr_set_si (res, 1, rnd_mode);
  if (mpfr_equal_p (x, __gmpfr_mone))
    return mpfr_set_si (res, (n & 1) == 0 ? 1 : -1, rnd_mode);

  /* P_0 = 1 */
  if (n == 0)
    return mpfr_set_si (res, 1, rnd_mode);

  /* P_1 = x */
  if (n == 1)
    {
      /* result is set to x. The ternary value of mpfr_set is returned */
      return mpfr_set (res, x, rnd_mode);
    }

  x_is_zero = MPFR_IS_ZERO (x);

  /* Pn(0) = 0 if n is odd */
  if (x_is_zero && (n & 1) == 1)
    {
      MPFR_SET_ZERO (res);

      /* Rationale: the a*x term of the Taylor expansion of P_n around x = 0
         (with n odd) has a > 0 for n mod 4 == 1, and a < 0 for n mod 4 == 3.
         Thus mpfr_legendre evaluates P_n(0) = +0 for n mod 4 == 1, and
         P_n(0) = -0 for n mod 4 == 3 */
      if ((n & 3) == 3)
        MPFR_SET_NEG (res);
      else
        MPFR_SET_POS (res);

      MPFR_RET (0);
    }

  res_prec = MPFR_PREC (res);

  /* Asymptotic expansion for small |x| and n >= 2.
     For the proof of the following bound, see algorithms.tex.
     Let t = n^2*x^2; the tail after the leading term is bounded by the
     following geometric series:
        |tail| <= |lead| * t / (1-t),
     where |lead| is c0 (for n even), or c1*x (for n odd) */
  if (!x_is_zero && n >= 2)
    {
      /* ex = MPFR_GET_EXP(x), such that 2^(ex-1) <= |x| < 2^ex;
         l2n = ceil(log2(n)), so n <= 2^l2n;
         thus, t = n^2*x^2 < (2^l2n)^2 * (2^ex)^2 = 2^{2*l2n+2*ex}.
         We define rho = 2*l2n+2*ex, therefore t < 2^{rho}.
         Note: we assume that rho is not going to overflow */
      ex = MPFR_GET_EXP (x);
      l2n = (mpfr_exp_t) MPFR_INT_CEIL_LOG2 (n);
      rho = 2 * l2n + 2 * ex;

      /* The bound err = -rho - 1 requires rho <= -2. In practice, require
         64 extra bits so the first rounding test usually succeeds */
      if (rho <= -2)
        {
          /* See algorithms.tex for the calculation of this error bound */
          err = -rho - 1;

          if (err >= (mpfr_exp_t) res_prec + MPFR_LEGENDRE_SMALL_X_GUARD)
            {
              inex = asymptotic_small_x (res, n, x, rnd_mode,
                                         &inex_round, err);

              /* if asymptotic_small_x sets inex_round to 0, then it cannot
                 round. In that case, our asymptotic expansion failed, so we
                 fall back to the usual Ziv loop. Otherwise, we return the
                 inex flag */
              if (inex_round)
                return inex;
            }
        }
    }

  /* Analyzing all the test cases where the result is not exact (inex != 0),
     we find that the average number of bits lost per iteration, i.e.,
     lost_bits/(n-1), is about 3.82. We thus add 4*n guard bits.
     For revision 94a5659, we have a total of 615575 such tests.
     With 4n+12 below, we get a probability of failure of 4.7%.
     With 4n+20, we get a probability of failure of 0.7%. */
  realprec = res_prec + 4 * n + 20;
  realprec += MPFR_INT_CEIL_LOG2 (realprec);

  MPFR_GROUP_INIT_5 (main, realprec,
                     p1, p2, pn, first_term, second_term);

  MPFR_ZIV_INIT (loop, realprec);
  for (;;)
    {
      /* p1 = x, p2 = 1 */
      inex = mpfr_set (p1, x, MPFR_RNDN);     /* p1 is a is algorithms.tex */
      /* If x is 0, set exp_p1 arbitrarily to 0 (necessarily n is even).
         In that case P_2n(0) = (-1)^n*(2n)!/(n!)^2/2^(2n). */
      mpfr_set_ui (p2, 1, MPFR_RNDN);         /* exact, p2 is b */

      /* In the loop:
           2^a_i is a bound on the absolute error on p1.
           2^b_i is a bound on the absolute error on p2.
         a_i and b_i come from the previous iteration, and initialized
         below for the first iteration (i = 2). */
      a_i = x_is_zero ? MPFR_EXP_MIN : MPFR_GET_EXP (p1) - realprec - 1;
      b_i = MPFR_EXP_MIN;

      for (i = 2; i <= n; i++)
        {
          MPFR_LOG_MSG (("i = %ld\n", i));
          MPFR_LOG_VAR (p1);
          MPFR_LOG_VAR (p2);

          /* FIXME: Extend the exponent range as usual? (There are currently
             no failures in the testsuite, but the tests may be incomplete.)
             Intermediate underflows are probably harmless, but must be
             taken into account in the error analysis. Final underflows
             (in mpfr_sub or mpfr_div_ui) might also be possible, and also
             need to be taken into account.
             Note: in the extended exponent range, underflows could occur
             only in very high precisions (thus only on 32-bit machines). */

          if (x_is_zero)
            {
              /* Special case when x = 0: first_term = 0, exact. */
              if ((i & 1) != 0)
                {
                  /* If i is odd, then p2 = 0, thus second_term = 0 too,
                     and pn = 0, exact. */
                  MPFR_ASSERTD (MPFR_IS_ZERO (p2));
                  MPFR_SET_ZERO (pn);
                  a_n = MPFR_EXP_MIN;
                  goto end_of_loop;
                }

              MPFR_SET_ZERO (first_term);
              h_i = MPFR_EXP_MIN;
            }
          else
            {
              mpfr_exp_t f_i;

              /* first_term = x * (2 * i - 1), with absolute error at step i
                 (denoted f_i in algorithms.tex)
                 bounded by f_i <= exp(first_term) - p - 1 */
              inex |= mpfr_mul_ui (first_term, x, 2 * i - 1, MPFR_RNDN);
              f_i = MPFR_GET_EXP (first_term) - realprec - 1;

              /* first_term = first_term * p1, with absolute error at step i
                 bounded by
                 h_i <= 2 + max(exp(first_term)-p-1, f_i+exp(p1),
                                MPFR_INT_CEIL_LOG2(2*i-1)+MPFR_GET_EXP(x)+a_i)
              */
              inex |= mpfr_mul (first_term, first_term, p1, MPFR_RNDN);
              h_i = 2 + MAX3 (MPFR_GET_EXP (first_term) - realprec - 1,
                              f_i + MPFR_GET_EXP (p1),
                              MPFR_INT_CEIL_LOG2 (2*i-1)
                              + MPFR_GET_EXP (x) + a_i);
            }

          log2_i_m1 = MPFR_INT_CEIL_LOG2 (i - 1);

          /* second_term = p2 * (i - 1), with absolute error at step i
             bounded by
             g_i <= max(exp(second_term)-p,
                        error(p2) + MPFR_INT_CEIL_LOG2(i-1)+1) */
          inex |= mpfr_mul_ui (second_term, p2, i - 1, MPFR_RNDN);
          g_i = MAX (MPFR_GET_EXP (second_term) - realprec,
                     b_i + log2_i_m1 + 1);

          /* pn = first_term - second_term, with absolute error at step i
             bounded by
             q_i <= 2 + max(exp(pn)-p-1,
                            error(first_term), error(second_term)) */
          inex |= mpfr_sub (pn, first_term, second_term, MPFR_RNDN);
          q_i = 2 + MAX3 (MPFR_GET_EXP (pn) - realprec - 1, h_i, g_i);

          /* pn = pn/i, with absolute error at step i
             bounded by a_i <= max(exp(pn)-p, q_i-MPFR_INT_CEIL_LOG2(i-1)+2) */
          inex |= mpfr_div_ui (pn, pn, i, MPFR_RNDN);
          a_n = MAX (MPFR_GET_EXP (pn) - realprec, q_i - log2_i_m1 + 2);

        end_of_loop:
          /* p2 = p1, p1 = pn */
          mpfr_swap (p2, p1); /* now p2 approximates P_{i-1}(x) */
          mpfr_swap (p1, pn); /* now p1 approximates P_i(x) */
          b_i = a_i;          /* 2^b_i is a bound on the absolute error on p2 */
          a_i = a_n;          /* 2^a_i is a bound on the absolute error on p1 */
        }

      /* Now p1 approximates P_n(x), and 2^a_i is a bound on its absolute error.
         Since ulp(p1) = 2^(EXP(p1)-realprec)
         we get the relative error is bounded by:
         2^(a_i - (EXP(p1) - realprec - 1)) */
      lost_bits = a_i - (MPFR_GET_EXP (p1) - realprec);

      /* if inex=0, then all the computation was exact, thus p1 is exactly
         P_n(x), otherwise we call MPFR_CAN_ROUND() to check if we can deduce
         the correct rounding */
      if (inex == 0 ||
          (lost_bits < realprec &&
           MPFR_CAN_ROUND (p1, realprec - lost_bits, res_prec, rnd_mode)))
        break;

      MPFR_ZIV_NEXT (loop, realprec);
      MPFR_GROUP_REPREC_5 (main, realprec,
                           p1, p2, pn, first_term, second_term);
    }
  MPFR_ZIV_FREE (loop);

  ternary_value = mpfr_set (res, p1, rnd_mode);

  MPFR_GROUP_CLEAR (main);

  return ternary_value;
}
