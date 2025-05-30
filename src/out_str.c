/* mpfr_out_str -- output a floating-point number to a stream

Copyright 1999, 2001-2002, 2004, 2006-2025 Free Software Foundation, Inc.
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

/* Warning! S should not contain "%". */
#define OUT_STR_RET(S)                          \
  do                                            \
    {                                           \
      int r;                                    \
      r = fprintf (stream, (S));                \
      return r < 0 ? 0 : r;                     \
    }                                           \
  while (0)

size_t
mpfr_out_str (FILE *stream, int base, size_t n_digits, mpfr_srcptr op,
              mpfr_rnd_t rnd_mode)
{
  char *s, *s0;
  size_t l;
  mpfr_exp_t e;
  int err, r;

  MPFR_ASSERTN ((base >= -36 && base <= -2) || (base >= 2 && base <= 62));

  if (MPFR_UNLIKELY (MPFR_IS_SINGULAR (op)))
    {
      if (MPFR_IS_NAN (op))
        OUT_STR_RET ("@NaN@");
      else if (MPFR_IS_INF (op))
        OUT_STR_RET (MPFR_IS_POS (op) ? "@Inf@" : "-@Inf@");
      else
        {
          MPFR_ASSERTD (MPFR_IS_ZERO (op));
          OUT_STR_RET (MPFR_IS_POS (op) ? "0" : "-0");
        }
    }

  s = mpfr_get_str (NULL, &e, base, n_digits, op, rnd_mode);

  s0 = s;
  /* for op=3.1416 we have s = "31416" and e = 1 */

  l = strlen (s) + 1; /* size of allocated block returned by mpfr_get_str
                         - may be incorrect, as only an upper bound? */

  /* outputs possible sign and significand */
  err = (*s == '-' && fputc (*s++, stream) == EOF)
    || fputc (*s++, stream) == EOF  /* leading digit */
    || fputc ((unsigned char) MPFR_DECIMAL_POINT, stream) == EOF
    || fputs (s, stream) == EOF;     /* trailing significand */
  mpfr_free_func (s0, l);
  if (MPFR_UNLIKELY (err))
    return 0;

  e--;  /* due to the leading digit */

  /* outputs exponent */
  r = fprintf (stream, (base <= 10 ?
                        "e%" MPFR_EXP_FSPEC "d" :
                        "@%" MPFR_EXP_FSPEC "d"), (mpfr_eexp_t) e);

  /* Check error from fprintf or integer overflow (wrapping) on size_t */
  if (MPFR_UNLIKELY (r < 0 || l + r < l))
    return 0;

  return l + r;
}
