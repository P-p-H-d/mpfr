/* mpfr_sub1 -- internal function to perform a "real" subtraction

Copyright 2001-2025 Free Software Foundation, Inc.
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

#include "mpfr-impl.h"

/* Warning! Because of possible aliases (e.g. from mpfr_fms), for the
   detection of reused arguments, do comparisons on the pointers to the
   significands instead of pointers to the MPFR numbers. */

/* compute sign(b) * (|b| - |c|), with |b| > |c|, diff_exp = EXP(b) - EXP(c)
   Returns 0 iff result is exact,
   a negative value when the result is less than the exact value,
   a positive value otherwise.
*/

int
mpfr_sub1 (mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, mpfr_rnd_t rnd_mode)
{
  int sign;
  mpfr_exp_t diff_exp, exp_a, exp_b;
  mpfr_prec_t cancel, cancel1;
  mp_size_t cancel2, an, bn, cn, cn0;
  mp_limb_t *ap, *bp, *cp;
  mp_limb_t carry, bb, cc;
  mpfr_prec_t aq, bq;
  int inexact, shift_b, shift_c, add_exp = 0;
  int cmp_low = 0; /* used for rounding to nearest: 0 if low(b) = low(c),
                      negative if low(b) < low(c), positive if low(b) > low(c) */
  int sh, k;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("b[%Pd]=%.*Rg c[%Pd]=%.*Rg rnd=%d",
      mpfr_get_prec (b), mpfr_log_prec, b,
      mpfr_get_prec (c), mpfr_log_prec, c, rnd_mode),
     ("a[%Pd]=%.*Rg", mpfr_get_prec (a), mpfr_log_prec, a));

  MPFR_TMP_MARK(marker);
  ap = MPFR_MANT(a);
  an = MPFR_LIMB_SIZE(a);

  (void) MPFR_GET_PREC (a);
  (void) MPFR_GET_PREC (b);
  (void) MPFR_GET_PREC (c);

  sign = mpfr_cmp2 (b, c, &cancel);

  if (MPFR_UNLIKELY(sign == 0))
    {
      MPFR_LOG_MSG (("sign=0\n", 0));
      if (rnd_mode == MPFR_RNDD)
        MPFR_SET_NEG (a);
      else
        MPFR_SET_POS (a);
      MPFR_SET_ZERO (a);
      MPFR_RET (0);
    }

  /* sign != 0, so that cancel has a valid value. */
  MPFR_LOG_MSG (("sign=%d cancel=%Pd\n", sign, cancel));
  MPFR_ASSERTD (cancel >= 0 && cancel <= MPFR_PREC_MAX);

  /*
   * If subtraction: sign(a) = sign * sign(b)
   * If addition: sign(a) = sign of the larger argument in absolute value.
   *
   * Both cases can be simplified in:
   * if (sign>0)
   *    if addition: sign(a) = sign * sign(b) = sign(b)
   *    if subtraction, b is greater, so sign(a) = sign(b)
   * else
   *    if subtraction, sign(a) = - sign(b)
   *    if addition, sign(a) = sign(c) (since c is greater)
   *      But if it is an addition, sign(b) and sign(c) are opposed!
   *      So sign(a) = - sign(b)
   */

  if (sign < 0) /* swap b and c so that |b| > |c| */
    {
      mpfr_srcptr t;
      MPFR_SET_OPPOSITE_SIGN (a,b);
      t = b; b = c; c = t;
    }
  else
    MPFR_SET_SAME_SIGN (a,b);

  if (MPFR_UNLIKELY (MPFR_IS_UBF (b) || MPFR_IS_UBF (c)))
    {
      exp_b = MPFR_UBF_GET_EXP (b);
      /* Early underflow detection. Rare, but a test is needed anyway
         since in the "MAX (aq, bq) + 2 <= diff_exp" branch, the exponent
         may decrease and MPFR_EXP_MIN would yield an integer overflow. */
      if (MPFR_UNLIKELY (exp_b < __gmpfr_emin - 1))
        {
          if (rnd_mode == MPFR_RNDN)
            rnd_mode = MPFR_RNDZ;
          return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
        }
      diff_exp = mpfr_ubf_diff_exp (b, c);
      MPFR_LOG_MSG (("UBF: exp_b=%" MPFR_EXP_FSPEC "d%s "
                     "diff_exp=%" MPFR_EXP_FSPEC "d%s\n",
                     (mpfr_eexp_t) exp_b,
                     exp_b == MPFR_EXP_MAX ? "=MPFR_EXP_MAX" : "",
                     (mpfr_eexp_t) diff_exp,
                     diff_exp == MPFR_EXP_MAX ? "=MPFR_EXP_MAX" : ""));
      /* If diff_exp == MPFR_EXP_MAX, the actual value can be larger,
         but anyway, since mpfr_exp_t >= mp_size_t, this will be the
         case c small below, and the exact value does not matter. */
      /* mpfr_set4 below used with MPFR_RNDF does not support UBF. */
      if (rnd_mode == MPFR_RNDF)
        rnd_mode = MPFR_RNDN;
    }
  else
    {
      exp_b = MPFR_GET_EXP (b);
      diff_exp = exp_b - MPFR_GET_EXP (c);
    }
  MPFR_ASSERTD (diff_exp >= 0);

  aq = MPFR_GET_PREC (a);
  bq = MPFR_GET_PREC (b);

  /* Check if c is too small.
     A more precise test is to replace 2 by
      (rnd == MPFR_RNDN) + mpfr_power2_raw (b)
      but it is more expensive and not very useful */
  if (MPFR_UNLIKELY (MAX (aq, bq) + 2 <= diff_exp))
    {
      MPFR_LOG_MSG (("case c small\n", 0));

      /* Remember, we can't have an exact result! */
      /*   A.AAAAAAAAAAAAAAAAA
         = B.BBBBBBBBBBBBBBB
          -                     C.CCCCCCCCCCCCC */
      /* A = S*ABS(B) +/- ulp(a) */

      /* since we can't have an exact result, for RNDF we can truncate b */
      if (rnd_mode == MPFR_RNDF)
        return mpfr_set4 (a, b, MPFR_RNDZ, MPFR_SIGN (a));

      exp_a = exp_b;  /* may be any out-of-range value due to UBF */
      MPFR_RNDRAW_EVEN (inexact, a, MPFR_MANT (b), bq,
                        rnd_mode, MPFR_SIGN (a),
                        if (exp_a != MPFR_EXP_MAX)
                          exp_a ++);
      MPFR_LOG_MSG (("inexact=%d\n", inexact));
      if (inexact == 0 &&
          /* a = b, but the exact value of b - c is a bit below. Then,
             except for directed rounding similar to toward zero and
             before overflow checking: a is the correctly rounded value
             and since |b| - |c| < |a|, the ternary value value is given
             by the sign of a. */
          ! MPFR_IS_LIKE_RNDZ (rnd_mode, MPFR_IS_NEG (a)))
        {
          MPFR_LOG_MSG (("c small, case 1\n", 0));
          inexact = MPFR_INT_SIGN (a);
        }
      else if (inexact != 0 &&
          /*   A.AAAAAAAAAAAAAA
             = B.BBBBBBBBBBBBBBB
              -                   C.CCCCCCCCCCCCC */
          /* It isn't exact, so PREC(b) > PREC(a) and the last
             PREC(b)-PREC(a) bits of b are not all zeros.
             Subtracting c from b will not have an effect on the rounding
             except in case of a midpoint in the round-to-nearest mode,
             when the even rounding was done away from zero instead of
             toward zero.
             In case of even rounding:
               1.BBBBBBBBBBBBBx10
             -                     1.CCCCCCCCCCCC
             = 1.BBBBBBBBBBBBBx01  Rounded to PREC(b)
             = 1.BBBBBBBBBBBBBx    Nearest / Rounded to PREC(a)
             Set gives:
               1.BBBBBBBBBBBBB0   if inexact == EVEN_INEX  (x == 0)
               1.BBBBBBBBBBBBB1+1 if inexact == -EVEN_INEX (x == 1)
             which means we get a wrong rounded result if x == 1,
             i.e. inexact == MPFR_EVEN_INEX (for positive numbers). */
               MPFR_LIKELY (inexact != MPFR_EVEN_INEX * MPFR_INT_SIGN (a)))
        {
          MPFR_LOG_MSG (("c small, case 2\n", 0));
          /* nothing to do */
        }
      else
        {
          /* We need to take the value preceding |a|. We can't use
             mpfr_nexttozero due to a possible out-of-range exponent.
             But this will allow us to have more specific code. */
          MPFR_LOG_MSG (("c small, case 3: correcting the value of a\n", 0));
          sh = (mpfr_prec_t) an * GMP_NUMB_BITS - aq;
          mpn_sub_1 (ap, ap, an, MPFR_LIMB_ONE << sh);
          if (MPFR_UNLIKELY (MPFR_LIMB_MSB (ap[an-1]) == 0))
            {
              exp_a --;
              /* The following is valid whether an = 1 or an > 1. */
              ap[an-1] |= MPFR_LIMB_HIGHBIT;
            }
          inexact = - MPFR_INT_SIGN (a);
        }
      /* The underflow case is possible only with UBF. The overflow case
         is also possible with normal FP due to rounding. */
      if (MPFR_UNLIKELY (exp_a > __gmpfr_emax))
        return mpfr_overflow (a, rnd_mode, MPFR_SIGN (a));
      if (MPFR_UNLIKELY (exp_a < __gmpfr_emin))
        {
          if (rnd_mode == MPFR_RNDN &&
              (exp_a < __gmpfr_emin - 1 ||
               (inexact * MPFR_INT_SIGN (a) >= 0 && mpfr_powerof2_raw (a))))
            rnd_mode = MPFR_RNDZ;
          return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
        }
      MPFR_SET_EXP (a, exp_a);
      MPFR_RET (inexact);
    }

  /* reserve a space to store b aligned with the result, i.e. shifted by
     (-cancel) % GMP_NUMB_BITS to the right */
  bn = MPFR_LIMB_SIZE (b);
  MPFR_UNSIGNED_MINUS_MODULO (shift_b, cancel);
  cancel1 = (cancel + shift_b) / GMP_NUMB_BITS;

  /* the high cancel1 limbs from b should not be taken into account */
  if (MPFR_UNLIKELY (shift_b == 0))
    {
      bp = MPFR_MANT(b); /* no need of an extra space */
      /* Ensure ap != bp */
      if (MPFR_UNLIKELY (ap == bp))
        {
          bp = MPFR_TMP_LIMBS_ALLOC (bn);
          MPN_COPY (bp, ap, bn);
        }
    }
  else
    {
      bp = MPFR_TMP_LIMBS_ALLOC (bn + 1);
      bp[0] = mpn_rshift (bp + 1, MPFR_MANT(b), bn++, shift_b);
    }

  /* reserve a space to store c aligned with the result, i.e. shifted by
     (diff_exp-cancel) % GMP_NUMB_BITS to the right */
  cn = MPFR_LIMB_SIZE (c);
  if (IS_POW2 (GMP_NUMB_BITS))
    shift_c = ((mpfr_uexp_t) diff_exp - cancel) % GMP_NUMB_BITS;
  else
    {
      /* The above operation does not work if diff_exp - cancel < 0. */
      shift_c = diff_exp - (cancel % GMP_NUMB_BITS);
      shift_c = (shift_c + GMP_NUMB_BITS) % GMP_NUMB_BITS;
    }
  MPFR_ASSERTD (shift_c >= 0 && shift_c < GMP_NUMB_BITS);

  if (MPFR_UNLIKELY(shift_c == 0))
    {
      cp = MPFR_MANT(c);
      /* Ensure ap != cp */
      if (ap == cp)
        {
          cp = MPFR_TMP_LIMBS_ALLOC (cn);
          MPN_COPY(cp, ap, cn);
        }
    }
  else
    {
      cp = MPFR_TMP_LIMBS_ALLOC (cn + 1);
      cp[0] = mpn_rshift (cp + 1, MPFR_MANT(c), cn++, shift_c);
    }

#if 0
  MPFR_LOG_MSG (("rnd=%s shift_b=%d shift_c=%d diffexp=%" MPFR_EXP_FSPEC
                 "d\n", mpfr_print_rnd_mode (rnd_mode), shift_b, shift_c,
                 (mpfr_eexp_t) diff_exp));
#endif

  MPFR_ASSERTD (ap != cp);
  MPFR_ASSERTD (bp != cp);

  /* here we have shift_c = (diff_exp - cancel) % GMP_NUMB_BITS,
        0 <= shift_c < GMP_NUMB_BITS
     thus we want cancel2 = ceil((cancel - diff_exp) / GMP_NUMB_BITS) */

  /* Possible optimization with a C99 compiler (i.e. well-defined
     integer division): if MPFR_PREC_MAX is reduced to
     ((mpfr_prec_t)((mpfr_uprec_t)(~(mpfr_uprec_t)0)>>1) - GMP_NUMB_BITS + 1)
     and diff_exp is of type mpfr_exp_t (no need for mpfr_uexp_t, since
     the sum or difference of 2 exponents must be representable, as used
     by the multiplication code), then the computation of cancel2 could
     be simplified to
       cancel2 = (cancel - (diff_exp - shift_c)) / GMP_NUMB_BITS;
     because cancel, diff_exp and shift_c are all non-negative and
     these variables are signed. */

  MPFR_ASSERTD (cancel >= 0);
  if (cancel >= diff_exp)
    /* Note that cancel is signed and will be converted to mpfr_uexp_t
       (type of diff_exp) in the expression below, so that this will
       work even if cancel is very large and diff_exp = 0. */
    cancel2 = (cancel - diff_exp + (GMP_NUMB_BITS - 1)) / GMP_NUMB_BITS;
  else
    cancel2 = - (mp_size_t) ((diff_exp - cancel) / GMP_NUMB_BITS);
  /* the high cancel2 limbs from b should not be taken into account */
#if 0
  MPFR_LOG_MSG (("cancel=%Pd cancel1=%Pd cancel2=%Pd\n",
                 cancel, cancel1, cancel2));
#endif

  /*               ap[an-1]        ap[0]
             <----------------+-----------|---->
             <----------PREC(a)----------><-sh->
 cancel1
 limbs        bp[bn-cancel1-1]
 <--...-----><----------------+-----------+----------->
  cancel2
  limbs       cp[cn-cancel2-1]                                    cancel2 >= 0
    <--...--><----------------+----------------+---------------->
                (-cancel2)                                        cancel2 < 0
                   limbs      <----------------+---------------->
  */

  /* first part: put in ap[0..an-1] the value of high(b) - high(c),
     where high(b) consists of the high an+cancel1 limbs of b,
     and high(c) consists of the high an+cancel2 limbs of c.
   */

  /* copy high(b) into a */
  if (MPFR_LIKELY(an + (mp_size_t) cancel1 <= bn))
    /* a: <----------------+-----------|---->
       b: <-----------------------------------------> */
      MPN_COPY (ap, bp + bn - (an + cancel1), an);
  else
    /* a: <----------------+-----------|---->
       b: <-------------------------> */
    if ((mp_size_t) cancel1 < bn) /* otherwise b does not overlap with a */
      {
        MPN_ZERO (ap, an + cancel1 - bn);
        MPN_COPY (ap + (an + cancel1 - bn), bp, bn - cancel1);
      }
    else
      MPN_ZERO (ap, an);

  /* subtract high(c) */
  if (MPFR_LIKELY(an + cancel2 > 0)) /* otherwise c does not overlap with a */
    {
      mp_limb_t *ap2;

      if (cancel2 >= 0)
        {
          if (an + cancel2 <= cn)
            /* a: <----------------------------->
               c: <-----------------------------------------> */
            mpn_sub_n (ap, ap, cp + cn - (an + cancel2), an);
          else
            /* a: <---------------------------->
               c: <-------------------------> */
            {
              ap2 = ap + an + (cancel2 - cn);
              if (cn > cancel2)
                mpn_sub_n (ap2, ap2, cp, cn - cancel2);
            }
        }
      else /* cancel2 < 0 */
        {
          mp_limb_t borrow;

          if (an + cancel2 <= cn)
            /* a: <----------------------------->
               c: <-----------------------------> */
            borrow = mpn_sub_n (ap, ap, cp + cn - (an + cancel2),
                                an + cancel2);
          else
            /* a: <---------------------------->
               c: <----------------> */
            {
              ap2 = ap + an + cancel2 - cn;
              borrow = mpn_sub_n (ap2, ap2, cp, cn);
            }
          ap2 = ap + an + cancel2;
          mpn_sub_1 (ap2, ap2, -cancel2, borrow);
        }
    }

  /* now perform rounding */
  sh = (mpfr_prec_t) an * GMP_NUMB_BITS - aq;
  /* last unused bits from a */
  carry = ap[0] & MPFR_LIMB_MASK (sh);
  ap[0] -= carry;

  if (rnd_mode == MPFR_RNDF)
    {
      inexact = 0;
      /* truncating is always correct since -1 ulp < low(b) - low(c) < 1 ulp */
      goto truncate;
    }
  else if (rnd_mode == MPFR_RNDN)
    {
      if (MPFR_LIKELY(sh))
        {
          /* can decide except when carry = 2^(sh-1) [middle]
             or carry = 0 [truncate, but cannot decide inexact flag] */
          if (carry > (MPFR_LIMB_ONE << (sh - 1)))
            goto add_one_ulp;
          else if ((0 < carry) && (carry < (MPFR_LIMB_ONE << (sh - 1))))
            {
              inexact = -1; /* result if smaller than exact value */
              goto truncate;
            }
          /* now carry = 2^(sh-1), in which case cmp_low=2,
             or carry = 0, in which case cmp_low=0 */
          cmp_low = (carry == 0) ? 0 : 2;
        }
    }
  else /* directed rounding: set rnd_mode to RNDZ iff toward zero */
    {
      if (MPFR_IS_RNDUTEST_OR_RNDDNOTTEST(rnd_mode, MPFR_IS_NEG(a)))
        rnd_mode = MPFR_RNDZ;

      if (carry)
        {
          if (rnd_mode == MPFR_RNDZ)
            {
              inexact = -1;
              goto truncate;
            }
          else /* round away */
            goto add_one_ulp;
        }
    }

  /* we have to consider the low (bn - (an+cancel1)) limbs from b,
     and the (cn - (an+cancel2)) limbs from c. */
  bn -= an + cancel1;
  cn0 = cn;
  cn -= an + cancel2;

#if 0
  MPFR_LOG_MSG (("last sh=%d bits from a are %Mu, bn=%Pd, cn=%Pd\n",
                 sh, carry, (mpfr_prec_t) bn, (mpfr_prec_t) cn));
#endif

  /* for rounding to nearest, we couldn't conclude up to here in the following
     cases:
     1. sh = 0, then cmp_low=0: we can either truncate, subtract one ulp
        or add one ulp: -1 ulp < low(b)-low(c) < 1 ulp
     2. sh > 0 but the low sh bits from high(b)-high(c) equal 2^(sh-1):
        -0.5 ulp <= -1/2^sh < low(b)-low(c)-0.5 < 1/2^sh <= 0.5 ulp
        we can't decide the rounding, in that case cmp_low=2:
        either we truncate and flag=-1, or we add one ulp and flag=1
     3. the low sh>0 bits from high(b)-high(c) equal 0: we know we have to
        truncate but we can't decide the ternary value, here cmp_low=0:
        -0.5 ulp <= -1/2^sh < low(b)-low(c) < 1/2^sh <= 0.5 ulp
        we always truncate and inexact can be any of -1,0,1
  */

  /* note: here cn might exceed cn0, in which case we consider a zero limb */
  for (k = 0; (bn > 0) || (cn > 0); k = 1)
    {
      /* if cmp_low < 0, we know low(b) - low(c) < 0
         if cmp_low > 0, we know low(b) - low(c) > 0
            (more precisely if cmp_low = 2, low(b) - low(c) = 0.5 ulp so far)
         if cmp_low = 0, so far low(b) - low(c) = 0 */

      /* get next limbs */
      bb = (bn > 0) ? bp[--bn] : 0;
      if ((cn > 0) && (cn-- <= cn0))
        cc = cp[cn];
      else
        cc = 0;

      /* cmp_low compares low(b) and low(c) */
      if (cmp_low == 0) /* case 1 or 3 */
        cmp_low = (bb < cc) ? -2+k : (bb > cc) ? 1 : 0;

      /* Case 1 for k=0 splits into 7 subcases:
         1a: bb > cc + half
         1b: bb = cc + half
         1c: 0 < bb - cc < half
         1d: bb = cc
         1e: -half < bb - cc < 0
         1f: bb - cc = -half
         1g: bb - cc < -half

         Case 2 splits into 3 subcases:
         2a: bb > cc
         2b: bb = cc
         2c: bb < cc

         Case 3 splits into 3 subcases:
         3a: bb > cc
         3b: bb = cc
         3c: bb < cc
      */

      /* the case rounding to nearest with sh=0 is special since one couldn't
         subtract above 1/2 ulp in the trailing limb of the result */
      if (rnd_mode == MPFR_RNDN && sh == 0 && k == 0) /* case 1 for k=0 */
        {
          mp_limb_t half = MPFR_LIMB_HIGHBIT;

          /* add one ulp if bb > cc + half
             truncate if cc - half < bb < cc + half
             sub one ulp if bb < cc - half
          */

          if (cmp_low < 0) /* bb < cc: -1 ulp < low(b) - low(c) < 0,
                              cases 1e, 1f and 1g */
            {
              if (cc >= half)
                cc -= half;
              else /* since bb < cc < half, bb+half < 2*half */
                bb += half;
              /* now we have bb < cc + half:
                 we have to subtract one ulp if bb < cc,
                 and truncate if bb > cc */
            }
          else if (cmp_low >= 0) /* bb >= cc, cases 1a to 1d */
            {
              if (cc < half)
                cc += half;
              else /* since bb >= cc >= half, bb - half >= 0 */
                bb -= half;
              /* now we have bb > cc - half: we have to add one ulp if bb > cc,
                 and truncate if bb < cc */
              if (cmp_low > 0)
                cmp_low = 2;
            }
        }

#if 0
      MPFR_LOG_MSG (("k=%d bb=%Mu cc=%Mu cmp_low=%d\n", k, bb, cc, cmp_low));
#endif

      if (cmp_low < 0) /* low(b) - low(c) < 0: either truncate or subtract
                          one ulp */
        {
          if (rnd_mode == MPFR_RNDZ)
            goto sub_one_ulp; /* set inexact=-1 */
          else if (rnd_mode != MPFR_RNDN) /* round away */
            {
              inexact = 1;
              goto truncate;
            }
          else /* round to nearest */
            {
              /* If cmp_low < 0 and bb > cc, then -0.5 ulp < low(b)-low(c) < 0,
                 whatever the value of sh.
                 If sh>0, then cmp_low < 0 implies that the initial neglected
                 sh bits were 0 (otherwise cmp_low=2 initially), thus the
                 weight of the new bits is less than 0.5 ulp too.
                 If k > 0 (and sh=0) this means that either the first neglected
                 limbs bb and cc were equal (thus cmp_low was 0 for k=0),
                 or we had bb - cc = -0.5 ulp or 0.5 ulp.
                 The last case is not possible here since we would have
                 cmp_low > 0 which is sticky.
                 In the first case (where we have cmp_low = -1), we truncate,
                 whereas in the 2nd case we have cmp_low = -2 and we subtract
                 one ulp.
              */
              if (bb > cc || sh > 0 || cmp_low == -1)
                {  /* -0.5 ulp < low(b)-low(c) < 0,
                      bb > cc corresponds to cases 1e and 1f1
                      sh > 0 corresponds to cases 3c and 3b3
                      cmp_low = -1 corresponds to case 1d3 (also 3b3) */
                  inexact = 1;
                  goto truncate;
                }
              else if (bb < cc) /* here sh = 0 and low(b)-low(c) < -0.5 ulp,
                                   this corresponds to cases 1g and 1f3 */
                goto sub_one_ulp;
              /* the only case where we can't conclude is sh=0 and bb=cc,
                 i.e., we have low(b) - low(c) = -0.5 ulp (up to now), thus
                 we don't know if we must truncate or subtract one ulp.
                 Note: for sh=0 we can't have low(b) - low(c) = -0.5 ulp up to
                 now, since low(b) - low(c) > 1/2^sh */
            }
        }
      else if (cmp_low > 0) /* 0 < low(b) - low(c): either truncate or
                               add one ulp */
        {
          if (rnd_mode == MPFR_RNDZ)
            {
              inexact = -1;
              goto truncate;
            }
          else if (rnd_mode != MPFR_RNDN) /* round away */
            goto add_one_ulp;
          else /* round to nearest */
            {
              if (bb > cc)
                {
                  /* if sh=0, then bb>cc means that low(b)-low(c) > 0.5 ulp,
                     and similarly when cmp_low=2 */
                  if (cmp_low == 2) /* cases 1a, 1b1, 2a and 2b1 */
                    goto add_one_ulp;
                  /* sh > 0 and cmp_low > 0: this implies that the sh initial
                     neglected bits were 0, and the remaining low(b)-low(c)>0,
                     but its weight is less than 0.5 ulp */
                  else /* 0 < low(b) - low(c) < 0.5 ulp, this corresponds to
                          cases 3a, 1d1 and 3b1 */
                    {
                      inexact = -1;
                      goto truncate;
                    }
                }
              else if (bb < cc) /* 0 < low(b) - low(c) < 0.5 ulp, cases 1c,
                                   1b3, 2b3 and 2c */
                {
                  inexact = -1;
                  goto truncate;
                }
              /* the only case where we can't conclude is bb=cc, i.e.,
                 low(b) - low(c) = 0.5 ulp (up to now), thus we don't know
                 if we must truncate or add one ulp. */
            }
        }
      /* after k=0, we cannot conclude in the following cases, we split them
         according to the values of bb and cc for k=1:
         1b. sh=0 and cmp_low = 1 and bb-cc = half [around 0.5 ulp]
             1b1. bb > cc: add one ulp, inex = 1
             1b2: bb = cc: cannot conclude
             1b3: bb < cc: truncate, inex = -1
         1d. sh=0 and cmp_low = 0 and bb-cc = 0 [around 0]
             1d1: bb > cc: truncate, inex = -1
             1d2: bb = cc: cannot conclude
             1d3: bb < cc: truncate, inex = +1
         1f. sh=0 and cmp_low = -1 and bb-cc = -half [around -0.5 ulp]
             1f1: bb > cc: truncate, inex = +1
             1f2: bb = cc: cannot conclude
             1f3: bb < cc: sub one ulp, inex = -1
         2b. sh > 0 and cmp_low = 2 and bb=cc [around 0.5 ulp]
             2b1. bb > cc: add one ulp, inex = 1
             2b2: bb = cc: cannot conclude
             2b3: bb < cc: truncate, inex = -1
         3b. sh > 0 and cmp_low = 0 [around 0]
             3b1. bb > cc: truncate, inex = -1
             3b2: bb = cc: cannot conclude
             3b3: bb < cc: truncate, inex = +1
      */
    }

  if ((rnd_mode == MPFR_RNDN) && cmp_low != 0)
    {
      /* even rounding rule */
      if ((ap[0] >> sh) & 1)
        {
          if (cmp_low < 0)
            goto sub_one_ulp;
          else
            goto add_one_ulp;
        }
      else
        inexact = (cmp_low > 0) ? -1 : 1;
    }
  else
    inexact = 0;
  goto truncate;

 sub_one_ulp: /* sub one unit in last place to a */
  mpn_sub_1 (ap, ap, an, MPFR_LIMB_ONE << sh);
  inexact = -1;
  goto end_of_sub;

 add_one_ulp: /* add one unit in last place to a */
  if (MPFR_UNLIKELY(mpn_add_1 (ap, ap, an, MPFR_LIMB_ONE << sh)))
    /* result is a power of 2: 11111111111111 + 1 = 1000000000000000 */
    {
      ap[an-1] = MPFR_LIMB_HIGHBIT;
      add_exp = 1;
    }
  inexact = 1; /* result larger than exact value */

 truncate:
  if (MPFR_UNLIKELY((ap[an-1] >> (GMP_NUMB_BITS - 1)) == 0))
    /* case 1 - epsilon */
    {
      ap[an-1] = MPFR_LIMB_HIGHBIT;
      add_exp = 1;
    }

 end_of_sub:
  /* we have to set MPFR_EXP(a) to MPFR_EXP(b) - cancel + add_exp, taking
     care of underflows/overflows in that computation, and of the allowed
     exponent range */
  MPFR_TMP_FREE (marker);
  if (MPFR_LIKELY(cancel))
    {
      cancel -= add_exp; /* OK: add_exp is an int equal to 0 or 1 */
      MPFR_ASSERTD (cancel >= 0);
      /* Detect an underflow case to avoid a possible integer overflow
         with UBF in the computation of exp_a. */
      if (MPFR_UNLIKELY (exp_b < __gmpfr_emin - 1))
        {
          if (rnd_mode == MPFR_RNDN)
            rnd_mode = MPFR_RNDZ;
          return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
        }
      exp_a = exp_b - cancel;
      /* The following assertion corresponds to a limitation of the MPFR
         implementation. It may fail with a 32-bit ABI and huge precisions,
         but this is practically impossible with a 64-bit ABI. This kind
         of issue is not specific to this function. */
      MPFR_ASSERTN (exp_b != MPFR_EXP_MAX || exp_a > __gmpfr_emax);
      if (MPFR_UNLIKELY (exp_a < __gmpfr_emin))
        {
        underflow:
          if (rnd_mode == MPFR_RNDN &&
              (exp_a < __gmpfr_emin - 1 ||
               (inexact >= 0 && mpfr_powerof2_raw (a))))
            rnd_mode = MPFR_RNDZ;
          return mpfr_underflow (a, rnd_mode, MPFR_SIGN(a));
        }
      /* We cannot have an overflow here, except for UBFs. Indeed:
         exp_a = exp_b - cancel + add_exp <= emax - 1 + 1 <= emax.
         For UBFs, we can have exp_b > emax. */
      if (exp_a > __gmpfr_emax)
        {
          MPFR_ASSERTD (exp_b > __gmpfr_emax);  /* since exp_b >= exp_a */
          return mpfr_overflow (a, rnd_mode, MPFR_SIGN (a));
        }
    }
  else /* cancel = 0: MPFR_EXP(a) <- MPFR_EXP(b) + add_exp */
    {
      /* in case cancel = 0, add_exp can still be 1, in case b is just
         below a power of two, c is very small, prec(a) < prec(b),
         and rnd=away or nearest */
      MPFR_ASSERTD (add_exp == 0 || add_exp == 1);
      /* Overflow iff exp_b + add_exp > __gmpfr_emax in Z, but we do
         a subtraction below to avoid a potential integer overflow in
         the case exp_b == MPFR_EXP_MAX. */
      if (MPFR_UNLIKELY (exp_b > __gmpfr_emax - add_exp))
        return mpfr_overflow (a, rnd_mode, MPFR_SIGN (a));
      exp_a = exp_b + add_exp;
      /* Warning: an underflow can happen for UBFs, for example when
         mpfr_add is called from mpfr_fmma or mpfr_fmms. */
      if (MPFR_UNLIKELY (exp_a < __gmpfr_emin))
        goto underflow;
      MPFR_ASSERTD (exp_a >= __gmpfr_emin);
    }
  MPFR_SET_EXP (a, exp_a);
  /* check that result is msb-normalized */
  MPFR_ASSERTD(ap[an-1] > ~ap[an-1]);
  MPFR_RET (inexact * MPFR_INT_SIGN (a));
}
