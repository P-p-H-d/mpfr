/* mpfr_ui_pow -- power of n function n^x

Copyright 2001-2004, 2006-2025 Free Software Foundation, Inc.
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

int
mpfr_ui_pow (mpfr_ptr y, unsigned long int n, mpfr_srcptr x, mpfr_rnd_t rnd_mode)
{
  mpfr_t t;
  int inexact;
  mp_limb_t tmp_mant[(sizeof (n) - 1) / MPFR_BYTES_PER_MP_LIMB + 1];
  MPFR_SAVE_EXPO_DECL (expo);

  MPFR_SAVE_EXPO_MARK (expo);
  MPFR_TMP_INIT1(tmp_mant, t, sizeof(n) * CHAR_BIT);
  inexact = mpfr_set_ui (t, n, MPFR_RNDN);
  MPFR_ASSERTD (inexact == 0);
  inexact = mpfr_pow (y, t, x, rnd_mode);
  MPFR_SAVE_EXPO_UPDATE_FLAGS (expo, __gmpfr_flags);
  MPFR_SAVE_EXPO_FREE (expo);
  return mpfr_check_range (y, inexact, rnd_mode);
}
