/*
 * frim.h -
 *
 * Definitions for FRiM (FRactal Iterative Method).
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of FRiM which is licensed under the MIT "Expat" License.
 *
 * Copyright (c) 2005-2020: Éric Thiébaut <https://github.com/emmt>
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _FRIM_H
#define _FRIM_H 1

#include <stdlib.h>

/*---------------------------------------------------------------------------*/
/* MACRO DEFINITIONS */

/*
 * _FRIM_EXTERN_C_BEGIN should be used at the beginning of your
 * declarations, so that C++ compilers don't mangle their names.  Use
 * _FRIM_EXTERN_C_END at the end of C declarations.
 */
#undef _FRIM_EXTERN_C_BEGIN
#undef _FRIM_EXTERN_C_END
#if defined(__cplusplus) || defined(c_plusplus)
#  if defined(__GNUC__) || defined(__clang__)
#    define restrict __restrict__
#  elif defined(_MSC_VER)
#    define restrict __restrict
#  else
#    define restrict
#  endif
#  define _FRIM_EXTERN_C_BEGIN extern "C" {
#  define _FRIM_EXTERN_C_END }
#else
#  define _FRIM_EXTERN_C_BEGIN
#  define _FRIM_EXTERN_C_END
#endif

/* A macro to allocate NBYTES bytes of dynamic memory. */
#define FRIM_ALLOC(nbytes)       malloc(nbytes)

/* A macro to free dynamic memory at address PTR. */
#define FRIM_FREE(ptr)           free(ptr)

/* A macro to allocate dynamic memory for NUMBER objects of type CLASS. */
#define FRIM_NEW(TYPE, NUMBER) \
                ((TYPE*)FRIM_ALLOC((NUMBER)*sizeof(TYPE)))

/* A macro to get the address offset of MEMBER in structure CLASS. */
#define FRIM_OFFSET_OF(TYPE, MEMBER) \
                ((char*)&((TYPE*)0)->MEMBER - (char*)0)

/* A macro to get the value at given offset from base address PTR. */
#define FRIM_FETCH_VALUE(TYPE, PTR, OFFSET) \
                (*(TYPE*)((char*)(PTR) + (OFFSET)))

/* A macro to compute smallest multiple of B greater or equal to A
   (arguments must be positive integers). */
#define FRIM_ROUND_UP(a, b)      ((((a) + (b) - 1)/(b))*(b))

/* A macro to get the largest of A and B. */
#define FRIM_MAX(a, b)           ((a) >= (b) ? (a) : (b))

/* A macro to get the smallest of A and B. */
#define FRIM_MIN(a, b)           ((a) <= (b) ? (a) : (b))

/* A macro to get the absolute value of A. */
#define FRIM_ABS(a)              ((a) >= 0 ? (a) : -(a))

/*
 * Utility macros: STRINGIFY takes an argument and wraps it in double quotation
 * marks (""), PASTE pastes two arguments.  Both perform macro expansion of
 * their arguments.
 */
#define  FRIM_VERBATIM(x)   x
#define  FRIM_STRINGIFY(x)  _FRIM_STRINGIFY(x)
#define _FRIM_STRINGIFY(x)  # x
#define  FRIM_PASTE(a,b)    _FRIM_PASTE_2(a, b)
#define _FRIM_PASTE_2(a,b)  a ## b

/*---------------------------------------------------------------------------*/

_FRIM_EXTERN_C_BEGIN

/*---------------------------------------------------------------------------*/
/* MESSAGES */

#define FRIM_SUCCESS        0  /* Integer code returned upon success. */
#define FRIM_FAILURE      (-1) /* Integer code returned upon failure. */
#define FRIM_MESSAGE_SIZE 128  /* Maximum size for error messages. */

/* Defining error messages along with their symbolic names in a macro in
   this way allows us to expand the macro in different contexts with
   confidence that the enumeration of symbolic names will map correctly
   onto the table of error messages.  */
#define _FRIM_ERROR_TABLE                                       \
    _FRIM_ERROR(NONE,          "OK - no error")                 \
    _FRIM_ERROR(NO_MEMORY,     "insufficient memory")           \
    _FRIM_ERROR(DIM_TOO_SMALL, "dimension too small")           \
    _FRIM_ERROR(BAD_DIM,       "bad dimension")                 \
    _FRIM_ERROR(BAD_NCOEFS,    "bad number of coefficients")    \
    _FRIM_ERROR(BAD_JOB,       "bad value for JOB parameter")   \
    _FRIM_ERROR(BAD_ARG,       "bad input argument")

/* Enumerate the symbolic error names. */
enum {
#define _FRIM_ERROR(ident, mesg) FRIM_PASTE(FRIM_ERROR_, ident),
    _FRIM_ERROR_TABLE
#undef _FRIM_ERROR
    FRIM_ERROR_MAX
};

/* Returns error corresponding to the given status. */
extern char const* frim_error_message(int const status);

/*---------------------------------------------------------------------------*/

/*
 * Various bit-flags (which can be combined with JOB).
 */

#define FRIM_JOB_DIRECT     0
#define FRIM_JOB_TRANSPOSE  1
#define FRIM_JOB_INVERSE    2
#define FRIM_SIX_VALUES     4 /* use 6 (instead of 4) input random values
                                 to generate the first 4 output random
                                 values */
#define FRIM_CLEAR          8 /* clear array (i.e. do not integrate
                                 gradient) */

/* Wavefront sensor model (No. 1). */
extern int wfs_model1_f(float*restrict dst, float const*restrict src,
                        size_t const dim, unsigned int flags, int* status);
extern int wfs_model1_d(double *restrict dst, double const*restrict src,
                        size_t const dim, unsigned int flags, int* status);

/* Lane et al. mid-point method for generating a Kolmogorov phase screen. */
extern int frim_lane_f(float*restrict dst, float const*restrict src,
                       size_t const dim, float const*restrict d, size_t nd,
                       unsigned int flags, int *status);
extern int frim_lane_d(double*restrict dst, double const*restrict src,
                       size_t const dim, double const*restrict d, size_t nd,
                       unsigned int flags, int* status);

/* 2D fractal generator for a given stationnary covariance. */
extern int frim_gen_2d_f(float*restrict dst, float const*restrict src,
                         size_t const dim, float const*restrict c,
                         size_t const nc, int job, int* status);
extern int frim_gen_2d_d(double *restrict dst, const double *restrict src,
                         size_t const dim, double const*restrict c,
                         size_t const nc, int job, int* status);

/* Bilinear fractal interpolator. */
extern int frim_int_2d_f(float*restrict dst, float const*restrict src,
                         size_t const dim, int job, int* status);
extern int frim_int_2d_d(double*restrict dst, double const*restrict src,
                         size_t const dim, int job, int* status);

/*---------------------------------------------------------------------------*/
_FRIM_EXTERN_C_END
#endif /* _FRIM_H */
