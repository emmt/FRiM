/*
 * frim.h --
 *
 *	Definitions for FRIM (FRactal Iterative Method).
 *
 *-----------------------------------------------------------------------------
 *
 *	Copyright (C) 2005-2006 Eric Thiébaut.
 *
 *	This file is part of FRIM (FRactal Iterative Method).
 *
 *	FRIM is free software; you can redistribute it and/or modify it
 *	under the terms of the GNU General Public License version 2 as
 *	published by the Free Software Foundation.
 *
 *	FRIM is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *	or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 *	License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with FRIM (file "COPYING" in the top source directory); if
 *	not, write to the Free Software Foundation, Inc., 51 Franklin St,
 *	Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *-----------------------------------------------------------------------------
 *
 *	$Id$
 *	$Log$
 *
 *-----------------------------------------------------------------------------
 */

#ifndef _FRIM_H
#define _FRIM_H 1

#include <stdlib.h>

/*---------------------------------------------------------------------------*/
/* MACRO DEFINITIONS */

/*
 * _FRIM_BEGIN_DECLS should be used at the beginning of your
 * declarations, so that C++ compilers don't mangle their names.  Use
 * _FRIM_END_DECLS at the end of C declarations.
 */
#undef _FRIM_BEGIN_DECLS
#undef _FRIM_END_DECLS
#if defined(__cplusplus) || defined(c_plusplus)
# define _FRIM_BEGIN_DECLS extern "C" {
# define _FRIM_END_DECLS }
#else
# define _FRIM_BEGIN_DECLS
# define _FRIM_END_DECLS
#endif

/* A macro to allocate NBYTES bytes of dynamic memory. */
#define FRIM_ALLOC(nbytes)       malloc(nbytes)

/* A macro to free dynamic memory at address PTR. */
#define FRIM_FREE(ptr)           free(ptr)

/* A macro to allocate dynamic memory for NUMBER objects of type CLASS. */
#define FRIM_NEW(class, number) \
		((class *)FRIM_ALLOC((number)*sizeof(class)))

/* A macro to get the address offset of MEMBER in structure CLASS. */
#define FRIM_OFFSET_OF(CLASS, MEMBER) \
		((char*)&((CLASS*)0)->MEMBER - (char*)0)

/* A macro to get the value at given offset from base address PTR. */
#define FRIM_FETCH_VALUE(TYPE, PTR, OFFSET) \
		(*(TYPE*)((char *)(PTR)+(OFFSET)))

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
 * Utility macros: STRINGIFY takes an argument and wraps it in "" (double
 * quotation marks), JOIN joins two arguments.  Both are capable of
 * performing macro expansion of their arguments.
 */
#define FRIM_VERBATIM(x) x
#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
# define FRIM_STRINGIFY(x)  FRIM_STRINGIFY1(x)
# define FRIM_STRINGIFY1(x) # x
# define FRIM_JOIN(a,b)     FRIM_JOIN1(a, b)
# define FRIM_JOIN1(a,b)    a ## b
#else
# define FRIM_STRINGIFY(x)  "x"
# define FRIM_JOIN(a,b)     FRIM_VERBATIM(a)/**/FRIM_VERBATIM(b)
#endif

/*---------------------------------------------------------------------------*/

_FRIM_BEGIN_DECLS

/*---------------------------------------------------------------------------*/
/* MESSAGES */

#define FRIM_SUCCESS        0  /* Integer code returned upon success. */
#define FRIM_FAILURE      (-1) /* Integer code returned upon failure. */
#define FRIM_MESSAGE_SIZE 128  /* Maximum size for error messages. */

/* Defining error messages along with their symbolic names in a macro in
   this way allows us to expand the macro in different contexts with
   confidence that the enumeration of symbolic names will map correctly
   onto the table of error messages.  */
#define _FRIM_ERROR_TABLE \
_FRIM_ERROR(NONE,          "OK - no error"), \
_FRIM_ERROR(NO_MEMORY,     "insufficient memory"), \
_FRIM_ERROR(DIM_TOO_SMALL, "dimension too small"), \
_FRIM_ERROR(BAD_DIM,       "bad dimension"), \
_FRIM_ERROR(BAD_NCOEFS,    "bad number of coefficients"), \
_FRIM_ERROR(BAD_JOB,       "bad value for JOB parameter"), \
_FRIM_ERROR(BAD_ARG,       "bad input argument")

/* Enumerate the symbolic error names. */
enum {
#define _FRIM_ERROR(ident, message) FRIM_JOIN(FRIM_ERROR_, ident)
  _FRIM_ERROR_TABLE,
#undef _FRIM_ERROR
  FRIM_ERROR_MAX
};

/* Returns error corresponding to the given status. */
extern const char *frim_error_message(const int status);

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
int wfs_model1_f(float dst[], const float src[], const size_t dim,
		 unsigned int flags, int *status);
int wfs_model1_d(double dst[], const double src[], const size_t dim,
		 unsigned int flags, int *status);

/* Lane et al. mid-point method for generating a Kolmogorov phase screen. */ 
extern int frim_lane_f(float dst[], const float src[], const size_t dim,
		       const float d[], size_t nd, unsigned int flags,
		       int *status);
extern int frim_lane_d(double dst[], const double src[], const size_t dim,
		       const double d[], size_t nd, unsigned int flags,
		       int *status);

/* 2D fractal generator for a given stationnary covariance. */
extern int frim_gen_2d_f(float dst[], const float src[], const size_t dim,
			 const float c[], const size_t nc, int job,
			 int *status);
extern int frim_gen_2d_d(double dst[], const double src[], const size_t dim,
			 const double c[], const size_t nc, int job,
			 int *status);

/* Bilinear fractal interpolator. */
extern int frim_int_2d_f(float dst[], const float src[], const size_t dim,
			 int job, int *status);
extern int frim_int_2d_d(double dst[], const double src[], const size_t dim,
			 int job, int *status);


/*---------------------------------------------------------------------------*/
_FRIM_END_DECLS
#endif /* _FRIM_H */
