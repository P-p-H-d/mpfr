/* __gmpfr_ceil_log2 - returns ceil(log(d)/log(2))

Copyright 1999-2004, 2006-2025 Free Software Foundation, Inc.
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

/* returns ceil(log(d)/log(2)) if d > 0,
   -1023 if d = +0,
   and floor(log(-d)/log(2))+1 if d < 0
*/
long
__gmpfr_ceil_log2 (double d)
{
  long exp;
#if _MPFR_IEEE_FLOATS
  union mpfr_ieee_double_extract x;

  x.d = d;
  /* The cast below is useless in theory, but let us not depend on the
     integer promotion rules (for instance, tcc is currently wrong). */
  exp = (long) x.s.exp - 1023;
  MPFR_ASSERTN (exp < 1023);  /* fail on infinities */
  x.s.exp = 1023; /* value for 1 <= d < 2 */
  if (x.d != 1.0) /* d: not a power of two? */
    exp++;
  return exp;
#else /* _MPFR_IEEE_FLOATS */
  double m;

  if (d < 0.0)
    return __gmpfr_floor_log2 (-d) + 1;
  else if (d == 0.0)
    return -1023;
  else if (d >= 1.0)
    {
      exp = 0;
      for (m = 1.0; m < d; m *= 2.0)
        exp++;
    }
  else
    {
      exp = 1;
      for (m = 1.0; m >= d; m *= 0.5)
        exp--;
    }
#endif /* _MPFR_IEEE_FLOATS */
  return exp;
}

