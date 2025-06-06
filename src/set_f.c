/* mpfr_set_f -- set a MPFR number from a GNU MPF number

Copyright 1999-2025 Free Software Foundation, Inc.
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
#define MPFR_NEED_MPF_INTERNALS
#include "mpfr-impl.h"

#ifndef MPFR_USE_MINI_GMP
int
mpfr_set_f (mpfr_ptr y, mpf_srcptr x, mpfr_rnd_t rnd_mode)
{
  mp_limb_t *my, *mx, *tmp;
  int cnt;
  mp_size_t sx, sy;
  int inexact, carry = 0;
  MPFR_TMP_DECL(marker);

  MPFR_LOG_FUNC
    (("rnd=%d", rnd_mode),
     ("y[%Pd]=%.*Rg inexact=%d", mpfr_get_prec (y), mpfr_log_prec, y,
      inexact));

  sx = ABSIZ(x); /* number of limbs of the mantissa of x */

  if (sx == 0) /* x is zero */
    {
      MPFR_SET_ZERO(y);
      MPFR_SET_POS(y);
      return 0; /* 0 is exact */
    }

  if (SIZ(x) * MPFR_FROM_SIGN_TO_INT(MPFR_SIGN(y)) < 0)
    MPFR_CHANGE_SIGN (y);

  sy = MPFR_LIMB_SIZE (y);
  my = MPFR_MANT(y);
  mx = PTR(x);

  count_leading_zeros(cnt, mx[sx - 1]);

  if (sy <= sx) /* we may have to round even when sy = sx */
    {
      mpfr_prec_t xprec = (mpfr_prec_t) sx * GMP_NUMB_BITS;

      MPFR_TMP_MARK(marker);
      tmp = MPFR_TMP_LIMBS_ALLOC (sx);
      if (cnt)
        mpn_lshift (tmp, mx, sx, cnt);
      else
        /* FIXME: we may avoid the copy here, and directly call mpfr_round_raw
           on mx instead of tmp */
        MPN_COPY (tmp, mx, sx);
      carry = mpfr_round_raw (my, tmp, xprec, (SIZ(x) < 0), MPFR_PREC(y),
                              rnd_mode, &inexact);
      if (MPFR_UNLIKELY(carry)) /* result is a power of two */
        my[sy - 1] = MPFR_LIMB_HIGHBIT;
      MPFR_TMP_FREE(marker);
    }
  else
    {
      if (cnt)
        mpn_lshift (my + sy - sx, mx, sx, cnt);
      else
        MPN_COPY (my + sy - sx, mx, sx);
      MPN_ZERO(my, sy - sx);
      /* no rounding necessary, since y has a larger mantissa */
      inexact = 0;
    }

  /* warning: EXP(x) * GMP_NUMB_BITS may exceed the maximal exponent */
  if (EXP(x) > 1 + (__gmpfr_emax - 1) / GMP_NUMB_BITS)
    {
      /* EXP(x) >= 2 + floor((__gmpfr_emax-1)/GMP_NUMB_BITS)
         EXP(x) >= 2 + (__gmpfr_emax - GMP_NUMB_BITS) / GMP_NUMB_BITS
                >= 1 + __gmpfr_emax / GMP_NUMB_BITS
         EXP(x) * GMP_NUMB_BITS >= __gmpfr_emax + GMP_NUMB_BITS
         Since 0 <= cnt <= GMP_NUMB_BITS-1, and 0 <= carry <= 1,
         we have then EXP(x) * GMP_NUMB_BITS - cnt + carry > __gmpfr_emax */
      return mpfr_overflow (y, rnd_mode, MPFR_SIGN (y));
    }
  else
    {
      /* Do not use MPFR_SET_EXP as the exponent may be out of range. */
      MPFR_EXP (y) = EXP (x) * GMP_NUMB_BITS - cnt + carry;
    }

  return mpfr_check_range (y, inexact, rnd_mode);
}
#endif
