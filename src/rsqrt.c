/* mpfr_rsqrt -- inverse square root following IEEE 754 (-Inf for -0)

Copyright 2008-2026 Free Software Foundation, Inc.
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

/* mpfr_rsqrt only differs from mpfr_rec_sqrt for -0,
   where mpfr_rec_sqrt returns +Inf (same result as for +0),
   and mpfr_rsqrt returns -Inf (following IEEE 754-2019) */
int
mpfr_rsqrt (mpfr_ptr r, mpfr_srcptr u, mpfr_rnd_t rnd_mode)
{
  MPFR_LOG_FUNC
    (("x[%Pd]=%.*Rg rnd=%d", mpfr_get_prec (u), mpfr_log_prec, u, rnd_mode),
     ("y[%Pd]=%.*Rg", mpfr_get_prec (r), mpfr_log_prec, r));

  /* IEEE 754-2019 says: rSqrt(-0) is -Inf and signals the divideByZero
     exception */
  if (MPFR_IS_ZERO(u) && MPFR_IS_NEG(u))
    {
      MPFR_SET_INF(r);
      MPFR_SET_NEG(r);
      MPFR_SET_DIVBY0 ();
      MPFR_RET(0); /* -Inf is exact */
    }

  /* for all other inputs, mpfr_rsqrt is identical to mpfr_rec_sqrt */
  return mpfr_rec_sqrt (r, u, rnd_mode);
}
