/* auxiliary functions for MPFR tests.

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

#ifndef __MPFR_TEST_H__
#define __MPFR_TEST_H__

/* Include config.h before using ANY configure macros if needed. */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

/* The no assertion request doesn't apply to the tests */
#if defined(MPFR_WANT_ASSERT)
# if MPFR_WANT_ASSERT < 0
#  undef MPFR_WANT_ASSERT
# endif
#endif

#include "mpfr-impl.h"

/* Avoid a GCC bug on Sparc, at least when using TLS. The MPFR library
 * itself is not affected, only a particular test. Normal code using
 * the MPFR library should not be affected either, as the bug occurs
 * when accessing __gmpfr_flags directly (and the public mpfr.h header
 * file does not provide any macro that accesses an internal variable
 * directly). So a workaround for the tests is the best solution.
 *
 * This bug, which could be observed under Debian with GCC 4.5.3 and
 * sparc-sun-solaris2.10 with GCC 5.5.0 when TLS and optimizations
 * are used[*], makes test programs using bad_cases() crash (SIGSEGV)
 * in this function at:
 *
 *   if (mpfr_nanflag_p () || mpfr_overflow_p () || mpfr_underflow_p ())
 *
 * Debugging shows that any attempt to access __gmpfr_flags directly
 * in the loop makes the program crash at this moment. This bug is not
 * present in the assembly code generated by -S, but is visible after a
 * normal compilation + link, when tracing the assembly code with gdb.
 * The workaround is to disable the macros from mpfr-impl.h that access
 * __gmpfr_flags directly. This bug may have been fixed in more recent
 * GCC versions, but it is safe to enable this workaround in all GCC
 * versions.
 *
 * [*] This is the default. Disabling TLS or recompiling the tests
 * without optimizations (-O0) makes the crash disappear.
 *
 * Mentions of these crashes:
 *   https://sympa.inria.fr/sympa/arc/mpfr/2011-10/msg00045.html [Debian]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2011-10/msg00055.html [Debian]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2011-12/msg00000.html [Solaris 10]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2011-12/msg00001.html [Solaris 10]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2011-12/msg00002.html
 *   https://sympa.inria.fr/sympa/arc/mpfr/2016-03/msg00061.html [Solaris 10]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2016-03/msg00063.html
 *   https://sympa.inria.fr/sympa/arc/mpfr/2020-06/msg00015.html [Solaris 10]
 *   https://sympa.inria.fr/sympa/arc/mpfr/2020-06/msg00020.html
 */
#if defined(__GNUC__) && defined (__sparc__)
# undef mpfr_underflow_p
# undef mpfr_overflow_p
# undef mpfr_nanflag_p
# undef mpfr_inexflag_p
# undef mpfr_erangeflag_p
# undef mpfr_divby0_p
#endif

#ifdef MPFR_TESTS_ABORT
# undef exit
# define exit(C) ((C) != 1 ? (exit)(C) : \
                  (fflush (stdout), fflush (stderr), abort ()))
#endif

#define STRINGIZE(S) #S
#define MAKE_STR(S) STRINGIZE(S)

/* In C (but not C++), mpfr_ptr and mpfr_srcptr arguments can be provided
   in a different pointer type, such as void *. For functions implemented
   as macros, the type conversion for the function parameters will not be
   done by the compiler, which means potential bugs in these implementations
   if we forget to take these unusual cases into account. So we need to test
   such arguments, in order to make sure that the arguments are converted to
   the expected type when needed.

   However, at least when the function is not implemented as a macro (which
   is the case when MPFR_USE_NO_MACRO is defined), such tests with void *
   arguments are not valid in C++; therefore, we will not do the cast to
   void * if the __cplusplus macro is defined. And with GCC 4.6+ compilers
   (and compatible), we will ignore the -Wc++-compat option around these
   tests.

   Note: in the future, inline functions could be used instead of macros,
   and such tests would become useless (except to detect compiler bugs).
*/
#if defined (__cplusplus)
#define VOIDP_CAST(X) (X)
#else
#define VOIDP_CAST(X) ((void *) (X))
/* Define IGNORE_CPP_COMPAT only for the GCC and Clang versions that
   support it.
   Note: GCC versions < 4.6 do not allow "#pragma GCC diagnostic" inside
   functions, and Clang on Windows (clang-cl) does not define __GNUC__.
   See https://sympa.inria.fr/sympa/arc/mpfr/2022-12/msg00007.html */
#if __MPFR_GNUC(4,6) || defined (__clang__)
#define IGNORE_CPP_COMPAT
#endif
#endif

#if defined (__cplusplus)
extern "C" {
#endif

/* generates a random long int, a random double,
   and corresponding seed initializing */
#define DBL_RAND() ((double) randlimb() / (double) MPFR_LIMB_MAX)

#define MINNORM 2.2250738585072013831e-308 /* 2^(-1022), smallest normalized */
#define MAXNORM 1.7976931348623157081e308 /* 2^(1023)*(2-2^(-52)) */

/* Generates a random rounding mode */
#define RND_RAND() ((mpfr_rnd_t) (randlimb() % MPFR_RND_MAX))

/* Ditto, excluding RNDF, assumed to be the last rounding mode */
#define RND_RAND_NO_RNDF() ((mpfr_rnd_t) (randlimb() % MPFR_RNDF))

/* Generates a random boolean (with type int, thanks to the "!= 0" test).
   Note: "& 1" is better than "% 2" for compilers with limited optimization,
   such as tcc 0.9.27. */
#define RAND_BOOL() ((randlimb() & 1) != 0)

/* Generates a random sign */
#define RAND_SIGN() (RAND_BOOL() ? MPFR_SIGN_POS : MPFR_SIGN_NEG)

/* Loop for all rounding modes */
#define RND_LOOP(_r) for((_r) = 0 ; (_r) < MPFR_RND_MAX ; (_r)++)

/* Loop for all rounding modes except RNDF (assumed to be the last one),
   which must be excluded from tests that rely on deterministic results. */
#define RND_LOOP_NO_RNDF(_r) for((_r) = 0 ; (_r) < MPFR_RNDF ; (_r)++)

/* Test whether two floating-point data have the same value,
   seen as an element of the set of the floating-point data
   (Level 2 in the IEEE 754-2008 standard). */
#define SAME_VAL(X,Y)                                                   \
  ((MPFR_IS_NAN (X) && MPFR_IS_NAN (Y)) ||                              \
   (mpfr_equal_p ((X), (Y)) && MPFR_INT_SIGN (X) == MPFR_INT_SIGN (Y)))

/* In the tests, mpfr_sgn was sometimes used incorrectly, for instance:
 *
 *   if (mpfr_cmp_ui (y, 0) || mpfr_sgn (y) < 0)
 *
 * to check that y is +0. This does not make sense since on 0, mpfr_sgn
 * yields 0, so that -0 would not be detected as an error. To make sure
 * that mpfr_sgn is not used incorrectly, we choose to fail when this
 * macro is used on a datum whose mathematical sign is not +1 or -1.
 * This feature is disabled when MPFR_TESTS_TSGN is defined, typically
 * in tsgn (to test mpfr_sgn itself).
 */
#ifndef MPFR_TESTS_TSGN
# undef mpfr_sgn
# define mpfr_sgn(x)                   \
  (MPFR_ASSERTN (! MPFR_IS_NAN (x)),   \
   MPFR_ASSERTN (! MPFR_IS_ZERO (x)),  \
   MPFR_SIGN (x))
#endif

#define FLIST mpfr_ptr, mpfr_srcptr, mpfr_rnd_t

/* For ofuf_thresholds in tests.c (argument threshold) */
#define POSOF (0+0)
#define POSUF (0+1)
#define NEGOF (2+0)
#define NEGUF (2+1)

int test_version (void);

/* Memory handling */
#define DEFAULT_MEMORY_LIMIT (1UL << 22)
extern size_t tests_memory_limit;
extern int tests_memory_disabled;
void *tests_allocate (size_t);
void *tests_reallocate (void *, size_t, size_t);
void tests_free (void *, size_t);
size_t tests_get_totalsize (void);
size_t tests_get_maxsize (void);
void tests_reset_maxsize (void);
void tests_memory_start (void);
void tests_memory_end (void);

void tests_start_mpfr (void);
void tests_end_mpfr (void);

void tests_expect_abort (void);
int tests_run_within_valgrind (void);

int mpfr_set_machine_rnd_mode (mpfr_rnd_t);
int have_subnorm_dbl (void);
int have_subnorm_flt (void);
void mpfr_test_init (void);
mp_limb_t randlimb (void);
unsigned long randulong (void);
long randlong (void);
void randseed (unsigned int);
void mpfr_random2 (mpfr_ptr, mp_size_t, mpfr_exp_t, gmp_randstate_t);
int ulp (double, double);
double dbl (double, int);
double Ulp (double);
int Isnan (double);
void d_trace (const char *, double);
void ld_trace (const char *, long double);
void n_trace (const char *, mp_limb_t *, mp_size_t);

FILE *src_fopen (const char *, const char *);
void set_emin (mpfr_exp_t);
void set_emax (mpfr_exp_t);
void tests_default_random (mpfr_ptr, int, mpfr_exp_t, mpfr_exp_t,
                           int);
void data_check (const char *, int (*) (FLIST), const char *);
void bad_cases (int (*)(FLIST), int (*)(FLIST),
                const char *, int, mpfr_exp_t, mpfr_exp_t,
                mpfr_prec_t, mpfr_prec_t, mpfr_prec_t, int);
void ofuf_thresholds (int (*)(FLIST), int (*)(FLIST), const char *,
                      mpfr_prec_t, mpfr_prec_t, int, unsigned int);
void flags_out (unsigned int);

int mpfr_cmp_str (mpfr_srcptr x, const char *, int, mpfr_rnd_t);
#define mpfr_cmp_str1(x,s) mpfr_cmp_str(x,s,10,MPFR_RNDN)
#define mpfr_set_str1(x,s) mpfr_set_str(x,s,10,MPFR_RNDN)

#define mpfr_cmp0(x,y) (MPFR_ASSERTN (!MPFR_IS_NAN (x) && !MPFR_IS_NAN (y)), mpfr_cmp (x,y))
#define mpfr_cmp_ui0(x,i) (MPFR_ASSERTN (!MPFR_IS_NAN (x)), mpfr_cmp_ui (x,i))
#define mpfr_cmp_si0(x,i) (MPFR_ASSERTN (!MPFR_IS_NAN (x)), mpfr_cmp_si (x,i))
#define mpfr_cmp_si_2exp0(x,i,e) (MPFR_ASSERTN (!MPFR_IS_NAN (x)), \
                                  mpfr_cmp_si_2exp (x,i,e))

/* define CHECK_EXTERNAL if you want to check mpfr against another library
   with correct rounding. You'll probably have to modify mpfr_print_raw()
   and/or test_add() below:
   * mpfr_print_raw() prints each number as "p m e" where p is the precision,
     m the mantissa (as a binary integer with sign), and e the exponent.
     The corresponding number is m*2^e. Example: "2 10 -6" represents
     2*2^(-6) with a precision of 2 bits.
   * test_add() outputs "b c a" on one line, for each addition a <- b + c.
     Currently it only prints such a line for rounding to nearest, when
     the inputs b and c are not NaN and/or Inf.
*/
#ifdef CHECK_EXTERNAL
static void
mpfr_print_raw (mpfr_srcptr x)
{
  printf ("%lu ", MPFR_PREC (x));
  if (MPFR_IS_NAN (x))
    {
      printf ("@NaN@");
      return;
    }

  if (MPFR_IS_NEG (x))
    printf ("-");

  if (MPFR_IS_INF (x))
    printf ("@Inf@");
  else if (MPFR_IS_ZERO (x))
    printf ("0 0");
  else
    {
      mp_limb_t *mx;
      mpfr_prec_t px;
      mp_size_t n;

      mx = MPFR_MANT (x);
      px = MPFR_PREC (x);

      for (n = (px - 1) / GMP_NUMB_BITS; ; n--)
        {
          mp_limb_t wd, t;

          MPFR_ASSERTN (n >= 0);
          wd = mx[n];
          for (t = MPFR_LIMB_HIGHBIT; t != 0; t >>= 1)
            {
              printf ((wd & t) == 0 ? "0" : "1");
              if (--px == 0)
                {
                  mpfr_exp_t ex;

                  ex = MPFR_GET_EXP (x);
                  MPFR_ASSERTN (ex >= LONG_MIN && ex <= LONG_MAX);
                  printf (" %ld", (long) ex - (long) MPFR_PREC (x));
                  return;
                }
            }
        }
    }
}
#endif

extern char *locale;
extern int tests_locale_enabled;

/* Random */
extern char             mpfr_rands_initialized;
extern gmp_randstate_t  mpfr_rands;

#undef RANDS
#define RANDS                                   \
  ((mpfr_rands_initialized ? 0                 \
    : (mpfr_rands_initialized = 1,             \
       gmp_randinit_default (mpfr_rands), 0)), \
   mpfr_rands)

#undef RANDS_CLEAR
#define RANDS_CLEAR()                   \
  do {                                  \
    if (mpfr_rands_initialized)        \
      {                                 \
        mpfr_rands_initialized = 0;    \
        gmp_randclear (mpfr_rands);    \
      }                                 \
  } while (0)

#if defined (__cplusplus)
}
#endif

/* With GCC, a macro "volatile" can be defined to test some special code
   in mpfr-impl.h (code for compilers that define such a macro), but the
   volatile keyword is necessary in some tests to avoid some GCC bugs.
   Thus we need to undef this macro (if defined). We do that at the end,
   so that mpfr-impl.h (included earlier) is not affected by this undef.
 */
#if defined(__GNUC__) && defined(volatile)
# undef volatile
#endif

#endif
