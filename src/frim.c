/*
 * frim.c --
 *
 *	FRIM (FRactal Iterative Method) main code.
 *
 *-----------------------------------------------------------------------------
 *
 *	Copyright (C) 2005-2006 Éric Thiébaut.
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

#ifndef _FRIM_CODE
#define _FRIM_CODE 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "frim.h"

#define RETURN_FAILURE(CODE) \
		if (status) *status = FRIM_JOIN(FRIM_ERROR_, CODE); \
		return FRIM_FAILURE
#define RETURN_SUCCESS \
		if (status) *status = FRIM_ERROR_NONE; \
		return FRIM_SUCCESS

#define real_t   float
#define suffix   _f
#include __FILE__

#define real_t   double
#define suffix   _d
#include __FILE__

#else /* _FRIM_CODE defined ---------------------------------------------*/

#define WFS_MODEL1   FRIM_JOIN(wfs_model1, suffix)
#define FRIM_GEN_2D  FRIM_JOIN(frim_gen_2d, suffix)
#define FRIM_INT_2D  FRIM_JOIN(frim_int_2d, suffix)
#define FRIM_LANE    FRIM_JOIN(frim_lane, suffix)

#ifdef WFS_MODEL1
int WFS_MODEL1(real_t dst[], const real_t src[], const size_t dim,
	       unsigned int flags, int *status)
{
  const real_t O_5 = 0.5;
  real_t wx, wy;
  size_t x, y, n;

  if (dim < 1) {
    RETURN_FAILURE(BAD_DIM);
  }
  n = dim - 1;

  if ((flags & FRIM_JOB_TRANSPOSE)) {
#define SRC(s,i1,i2) src[((i2)*n + (i1))*2 + (s)]
#define DST(i1,i2)   dst[(i2)*dim + (i1)]
    if ((flags & FRIM_CLEAR)) {
      memset(dst, 0, dim*dim*sizeof(*dst));
    }
    for (y=0 ; y<n ; ++y) {
      for (x=0 ; x<n ; ++x) {
	wx = SRC(0,x,y)*O_5;
	wy = SRC(1,x,y)*O_5;
	DST( x , y ) -= wx + wy;
	DST(x+1, y ) += wx - wy;
	DST( x ,y+1) -= wx - wy;
	DST(x+1,y+1) += wx + wy;
      }
    }
#undef SRC
#undef DST
  } else {
    /* DST is a 2×(DIM-1)×(DIM-1) array */
#define SRC(i1,i2)   src[(i2)*dim + (i1)]
#define DST(s,i1,i2) dst[((i2)*n + (i1))*2 + (s)]
    for (y=0 ; y<n ; ++y) {
      for (x=0 ; x<n ; ++x) {
	DST(0,x,y) = (SRC(x+1,y+1) + SRC(x+1,y) - SRC(x,y+1) - SRC(x,y))*O_5;
	DST(1,x,y) = (SRC(x+1,y+1) - SRC(x+1,y) + SRC(x,y+1) - SRC(x,y))*O_5;
      }
    }
#undef SRC
#undef DST
  }
  RETURN_SUCCESS;
}
#endif /* WFS_MODEL1 */

/*
 *  1--+--+--+--1  1--+--+--+--1  1--+--3--+--1
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  +--+--+--+--+  +--+--2--+--+  3--+--2--+--3
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  1--+--+--+--1  1--+--+--+--1  1--+--3--+--1
 *
 *  1--+--3--+--1  1--+--3--+--1  1--6--3--6--1
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  +--4--+--4--+  +--4--5--4--+  6--4--5--4--6
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  3--+--2--+--3  3--5--2--5--3  3--5--2--5--3
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  +--4--+--4--+  +--4--5--4--+  6--4--5--4--6
 *  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 *  1--+--3--+--1  1--+--3--+--1  1--6--3--6--1
 */

/*
 * Assuming STP is the step size between adjacent points, the array C
 * is:
 *   C[0] = Var = Cov(0)
 *   C[1] = Cov(STP)
 *   C[2] = Cov(STP*sqrt(2))
 *   ...
 *   C[k] = Cov(STP*2^((k-1)/2))  for k = 1,..., NC - 1
 *   ...
 *   C[NC - 2] = Cov(STP*(DIM - 1))
 *   C[NC - 1] = Cov(STP*sqrt(2)*(DIM - 1))
 *
 * NC - Number of elements for array of precomputed covariances.
 *
 * DIM = 2^P + 1
 * NC = 2*P + 3
 */

#ifdef FRIM_GEN_2D
int FRIM_GEN_2D(real_t dst[], const real_t src[], const size_t dim,
		const real_t c[], const size_t nc, int job, int *status)
{
  real_t c0, c1, c2, c3;
  real_t w1, w2, w3, w4, w5, q;
  real_t p1, p2, p3, p4;
  size_t n, p, k, x, y, s, h;

  /* Check the arguments. */
  n = dim - 1;
  for (p = 0 ; 1<<p != n ; ) {
    if (++p >= (8*sizeof(size_t))) {
      RETURN_FAILURE(BAD_DIM);
    }
  }
  if (nc != 2*p + 3) {
      RETURN_FAILURE(BAD_NCOEFS);
  }

  /* Copy SRC in DST array to perform 'in-place' operation. */
  if (src && src != dst) {
    memcpy(dst, src, dim*dim*sizeof(src[0]));
  }

#define A(i1,i2) dst[(i2)*dim + (i1)]

  n = dim - 1;
  c0 = c[0]; /* C0 <- Cov(0) = Var */

  if (job == 0) {

    /*******************************
     **                           **
     **  APPLY FORWARD TRANSFORM  **
     **                           **
     *******************************/

    /* Walk from largest scales to smallest ones. */
    k = nc;

    /* Outermost 4 corners. */
    c2 = c[--k]; /* C2 <- Cov(sqrt(2)*(DIM - 1)) */
    c1 = c[--k]; /* C1 <- Cov(DIM - 1) */
    w1 = sqrt(c0 + 2.0*c1 + c2)/2.0;
    w2 = sqrt(c0 - 2.0*c1 + c2)/2.0;
    w3 = sqrt((c0 - c2)/2.0);
    p1 = w1*A(0,0);
    p2 = w2*A(n,0);
    p3 = w3*A(0,n);
    p4 = w3*A(n,n);
    A(0,0) = p1 - p2 - p3;
    A(n,0) = p1 + p2 - p4;
    A(0,n) = p1 + p2 + p4;
    A(n,n) = p1 - p2 + p3;

    for (s=n ; s>=2 ; s/=2) {
      /* Shift scale factors. */
      h = s/2;
      c3 = c2;
      c2 = c1;
      c1 = c[--k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  A(x,y) = (A(x-h,y-h) +
		    A(x+h,y-h) +
		    A(x-h,y+h) +
		    A(x+h,y+h))*w1 + A(x,y)*w2;
	}
      }

      /* Shift scale factors. */
      c3 = c2;
      c2 = c1;
      c1 = c[--k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      q = c1/((c0 + c3)*c0 - 2.0*c2*c2);
      w3 = (c0 - c2)*q;
      w4 = (c0 - 2.0*c2 + c3)*q;
      w5 = sqrt(c0 - (w3 + w4 + w3)*c1);

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	A(x,0) = (A(x-h,0) + A(x+h,0))*w3 + A(x,h)*w4 + A(x,0)*w5;
      }
      for (y=h ; y<n ; y+=h) {
	A(0,y) = (A(0,y-h) + A(0,y+h))*w3 + A(h,y)*w4 + A(0,y)*w5;
	for (x=s ; x<n ; x+=s) {
	  A(x,y) = (A( x ,y-h) +
		    A(x-h, y ) +
		    A(x+h, y ) +
		    A( x ,y+h))*w1 + A(x,y)*w2;
	}
	A(n,y) = (A(n,y-h) + A(n,y+h))*w3 + A(n-h,y)*w4 + A(n,y)*w5;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  A(x,y) = (A( x ,y-h) +
		    A(x-h, y ) +
		    A(x+h, y ) +
		    A( x ,y+h))*w1 + A(x,y)*w2;
	}
      }
      for (x=h ; x<n ; x+=s) {
	A(x,n) = A(x,n-h)*w4 + (A(x-h,n) + A(x+h,n))*w3 + A(x,n)*w5;
      }
    }

  } else if (job == 1) {

    /*******************************
     **                           **
     **  APPLY TRANSPOSE OPERATOR **
     **                           **
     *******************************/

    /* Walk from smallest scales to largest ones. */
    k = 0;
    c2 = c[++k];
    c3 = c[++k];

    for (s=2 ; s<=n ; s+=s) {

      /* Half-step size. */
      h = s/2;

      /* Shift scale factors. */
      c1 = c2;
      c2 = c3;
      c3 = c[++k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      q = c1/((c0 + c3)*c0 - 2.0*c2*c2);
      w3 = (c0 - c2)*q;
      w4 = (c0 - 2.0*c2 + c3)*q;
      w5 = sqrt(c0 - (w3 + w4 + w3)*c1);

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	q = A(x,0);
	A( x ,0)  = w5*q;
	A(x-h,0) += w3*q;
	A(x+h,0) += w3*q;
	A( x ,h) += w4*q;
      }
      for (y=h ; y<n ; y+=h) {
	q = A(0,y);
	A(0,y)    = w5*q;
	A(0,y-h) += w3*q;
	A(0,y+h) += w3*q;
	A(h,y)   += w4*q;
	for (x=s ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x , y )  = w2*q;
	  A( x ,y-h) += w1*q;
	  A(x-h, y ) += w1*q;
	  A(x+h, y ) += w1*q;
	  A( x ,y+h) += w1*q;
	}
	q = A(n,y);
	A( n , y )  = w5*q;
	A( n ,y-h) += w3*q;
	A( n ,y+h) += w3*q;
	A(n-h, y ) += w4*q;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x,  y )  = w2*q;
	  A( x ,y-h) += w1*q;
	  A(x-h, y ) += w1*q;
	  A(x+h, y ) += w1*q;
	  A( x ,y+h) += w1*q;
	}
      }
      for (x=h ; x<n ; x+=s) {
	q = A(x,n);
	A( x , n )  = w5*q;
	A( x ,n-h) += w4*q;
	A(x-h, n ) += w3*q;
	A(x+h, n ) += w3*q;
      }

      /* Shift scale factors. */
      c1 = c2;
      c2 = c3;
      c3 = c[++k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x , y )  = w2*q;
	  A(x-h,y-h) += w1*q;
	  A(x+h,y-h) += w1*q;
	  A(x-h,y+h) += w1*q;
	  A(x+h,y+h) += w1*q;
	}
      }
    }

    /* Outermost 4 corners. */
    c1 = c2; /* C1 <- Cov(DIM - 1) */
    c2 = c3; /* C2 <- Cov(sqrt(2)*(DIM - 1)) */
    w1 = sqrt(c0 + 2.0*c1 + c2)*0.5;
    w2 = sqrt(c0 - 2.0*c1 + c2)*0.5;
    w3 = sqrt((c0 - c2)*0.5);
    p1 = A(0,0);
    p2 = A(n,0);
    p3 = A(0,n);
    p4 = A(n,n);
    A(0,0) = w1*(p1 + p2 + p3 + p4);
    A(n,0) = w2*(p2 - p1 + p3 - p4);
    A(0,n) = w3*(p4 - p1);
    A(n,n) = w3*(p3 - p2);

  } else if (job == 2) {

    /*******************************
     **                           **
     **  APPLY INVERSE TRANSFORM  **
     **                           **
     *******************************/

    /* Walk from smallest scales to largest ones. */
    k = 0;
    c2 = c[++k];
    c3 = c[++k];

    for (s=2 ; s<=n ; s+=s) {

      /* Half-step size. */
      h = s/2;

      /* Shift scale factors. */
      c1 = c2;
      c2 = c3;
      c3 = c[++k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      w1 /= w2;
      w2 = 1.0/w2;
      q = c1/((c0 + c3)*c0 - 2.0*c2*c2);
      w3 = (c0 - c2)*q;
      w4 = (c0 - 2.0*c2 + c3)*q;
      w5 = sqrt(c0 - (w3 + w4 + w3)*c1);
      w3 /= w5;
      w4 /= w5;
      w5 = 1.0/w5;

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	A(x,0) = A(x,0)*w5 - (A(x-h,0) + A(x+h,0))*w3 - A(x,h)*w4;
      }
      for (y=h ; y<n ; y+=h) {
	A(0,y) = A(0,y)*w5 - (A(0,y-h) + A(0,y+h))*w3 - A(h,y)*w4;
	for (x=s ; x<n ; x+=s) {
	  A(x,y) = A(x,y)*w2 - (A( x ,y-h) +
				A(x-h, y ) +
				A(x+h, y ) +
				A( x ,y+h))*w1;
	}
	A(n,y) = A(n,y)*w5 - (A(n,y-h) + A(n,y+h))*w3 - A(n-h,y)*w4;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  A(x,y) = A(x,y)*w2 - (A( x ,y-h) +
				A(x-h, y ) +
				A(x+h, y ) +
				A( x ,y+h))*w1;
	}
      }
      for (x=h ; x<n ; x+=s) {
	A(x,n) = A(x,n)*w5 - A(x,n-h)*w4 - (A(x-h,n) + A(x+h,n))*w3;
      }

      /* Shift scale factors. */
      c1 = c2;
      c2 = c3;
      c3 = c[++k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      w1 /= w2;
      w2 = 1.0/w2;

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  A(x,y) = A(x,y)*w2 - (A(x-h,y-h) +
				A(x+h,y-h) +
				A(x-h,y+h) +
				A(x+h,y+h))*w1;
	}
      }
    }

    /* Outermost 4 corners. */
    c1 = c2; /* C1 <- Cov(DIM - 1) */
    c2 = c3; /* C2 <- Cov(sqrt(2)*(DIM - 1)) */
#if 0 /* unused because seems to produce larger rounding errors */
    w1 = sqrt(0.25/(c0 + 2.0*c1 + c2));
    w2 = sqrt(0.25/(c0 - 2.0*c1 + c2));
    w3 = sqrt(0.5/(c0 - c2));
#else
    w1 = 1.0/sqrt((c0 + 2.0*c1 + c2)*4.0);
    w2 = 1.0/sqrt((c0 - 2.0*c1 + c2)*4.0);
    w3 = 1.0/sqrt((c0 - c2)*2.0);
#endif
    p1 = A(0,0);
    p2 = A(n,0);
    p3 = A(0,n);
    p4 = A(n,n);
    A(0,0) = (p1 + p2 + p3 + p4)*w1;
    A(n,0) = (p2 - p1 + p3 - p4)*w2;
    A(0,n) = (p4 - p1)*w3;
    A(n,n) = (p3 - p2)*w3;

  } else if (job == 3) {

    /*****************************************
     **                                     **
     **  APPLY INVERSE TRANSPOSE TRANSFORM  **
     **                                     **
     *****************************************/

    /* Walk from largest scales to smallest ones. */
    k = nc;

    /* Outermost 4 corners. */
    c2 = c[--k]; /* C2 <- Cov(sqrt(2)*(DIM - 1)) */
    c1 = c[--k]; /* C1 <- Cov(DIM - 1) */
#if 0
    w1 = sqrt((c0 + 2.0*c1 + c2)*4.0);
    w2 = sqrt((c0 - 2.0*c1 + c2)*4.0);
#endif
    p1 = A(0,0)/sqrt((c0 + 2.0*c1 + c2)*4.0);
    p2 = A(n,0)/sqrt((c0 - 2.0*c1 + c2)*4.0);
    w3 = sqrt((c0 - c2)*2.0);
    p3 = A(0,n)/w3;
    p4 = A(n,n)/w3;
    A(0,0) = p1 - p2 - p3;
    A(n,0) = p1 + p2 - p4;
    A(0,n) = p1 + p2 + p4;
    A(n,n) = p1 - p2 + p3;

    for (s=n ; s>=2 ; s/=2) {
      /* Shift scale factors. */
      h = s/2;
      c3 = c2;
      c2 = c1;
      c1 = c[--k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      w1 /= w2;
      w2 = 1.0/w2;

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x , y )  = w2*q;
	  A(x-h,y-h) -= w1*q;
	  A(x+h,y-h) -= w1*q;
	  A(x-h,y+h) -= w1*q;
	  A(x+h,y+h) -= w1*q;
	}
      }

      /* Shift scale factors. */
      c3 = c2;
      c2 = c1;
      c1 = c[--k];
      w1 = c1/(c0 + 2.0*c2 + c3);
      w2 = sqrt(c0 - 4.0*c1*w1);
      w1 /= w2;
      w2 = 1.0/w2;
      q = c1/((c0 + c3)*c0 - 2.0*c2*c2);
      w3 = (c0 - c2)*q;
      w4 = (c0 - 2.0*c2 + c3)*q;
      w5 = sqrt(c0 - (w3 + w4 + w3)*c1);
      w3 /= w5;
      w4 /= w5;
      w5 = 1.0/w5;

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	q = A(x,0);
	A( x ,0)  = w5*q;
	A(x-h,0) -= w3*q;
	A(x+h,0) -= w3*q;
	A( x ,h) -= w4*q;
      }
      for (y=h ; y<n ; y+=h) {
	q = A(0,y);
	A(0,y) = w5*q;
	A(0,y-h) -= w3*q;
	A(0,y+h) -= w3*q;
	A(h,y)   -= w4*q;
	for (x=s ; x<n ; x+=s) {
	  q = A(x,y);
	  A(x,y) = w2*q;
	  A( x ,y-h) -= w1*q;
	  A(x-h, y ) -= w1*q;
	  A(x+h, y ) -= w1*q;
	  A( x ,y+h) -= w1*q;
	}
	q = A(n,y);
	A(n,y) = w5*q;
	A(n,y-h) -= w3*q;
	A(n,y+h) -= w3*q;
	A(n-h,y) -= w4*q;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A(x,y) = w2*q;
	  A( x ,y-h) -= w1*q;
	  A(x-h, y ) -= w1*q;
	  A(x+h, y ) -= w1*q;
	  A( x ,y+h) -= w1*q;
	}
      }
      for (x=h ; x<n ; x+=s) {
	q = A(x,n);
	A(x,n) = w5*q;
	A(x,n-h) -= w4*q;
	A(x-h,n) -= w3*q;
	A(x+h,n) -= w3*q;
      }
    }

  } else {

    RETURN_FAILURE(BAD_JOB);

  }

#undef A

  RETURN_SUCCESS;
}
#endif /* FRIM_GEN_2D */

#ifdef FRIM_INT_2D
int FRIM_INT_2D(real_t dst[], const real_t src[], const size_t dim, int job,
		int *status)
{
  const real_t p25 = 1.0/4.0;
  const real_t p33 = 1.0/3.0;
  real_t q;
  size_t n, p, x, y, s, h;

  /* Check the arguments. */
  n = dim - 1;
  for (p = 0 ; 1<<p != n ; ) {
    if (++p >= (8*sizeof(size_t))) {
      RETURN_FAILURE(BAD_DIM);
    }
  }

  /* Copy SRC in DST array to perform 'in-place' operation. */
  if (src && src != dst) {
    memcpy(dst, src, dim*dim*sizeof(src[0]));
  }

#define A(i1,i2) dst[(i2)*dim + (i1)]

  if (job == 0) {

    /*******************************
     **                           **
     **  APPLY FORWARD TRANSFORM  **
     **                           **
     *******************************/

    /* Walk from largest scales to smallest ones.
       Outermost 4 corners are left unchanged. */
    for (s=n ; s>=2 ; s/=2) {

      /* Half-step size. */
      h = s/2;

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  A(x,y) += (A(x-h,y-h) + A(x+h,y-h) + A(x-h,y+h) + A(x+h,y+h))*p25;
	}
      }

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	A(x,0) += (A(x-h,0) + A(x+h,0) + A(x,h))*p33;
      }
      for (y=h ; y<n ; y+=h) {
	A(0,y) += (A(0,y-h) + A(0,y+h) + A(h,y))*p33;
	for (x=s ; x<n ; x+=s) {
	  A(x,y) += (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
	}
	A(n,y) += (A(n,y-h) + A(n,y+h) + A(n-h,y))*p33;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  A(x,y) += (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
	}
      }
      for (x=h ; x<n ; x+=s) {
	A(x,n) += (A(x,n-h) + A(x-h,n) + A(x+h,n))*p33;
      }
    }

  } else if (job == 1) {

    /*******************************
     **                           **
     **  APPLY TRANSPOSE OPERATOR **
     **                           **
     *******************************/

    /* Walk from smallest scales to largest ones. */
    for (s=2 ; s<=n ; s+=s) {

      /* Half-step size. */
      h = s/2;

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	q = A(x,0);
	A(x-h,0) += p33*q;
	A(x+h,0) += p33*q;
	A( x ,h) += p33*q;
      }
      for (y=h ; y<n ; y+=h) {
	q = A(0,y);
	A(0,y-h) += p33*q;
	A(0,y+h) += p33*q;
	A(h,y)   += p33*q;
	for (x=s ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x ,y-h) += p25*q;
	  A(x-h, y ) += p25*q;
	  A(x+h, y ) += p25*q;
	  A( x ,y+h) += p25*q;
	}
	q = A(n,y);
	A( n ,y-h) += p33*q;
	A( n ,y+h) += p33*q;
	A(n-h, y ) += p33*q;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A( x ,y-h) += p25*q;
	  A(x-h, y ) += p25*q;
	  A(x+h, y ) += p25*q;
	  A( x ,y+h) += p25*q;
	}
      }
      for (x=h ; x<n ; x+=s) {
	q = A(x,n);
	A( x ,n-h) += p33*q;
	A(x-h, n ) += p33*q;
	A(x+h, n ) += p33*q;
      }

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y);
	  A(x-h,y-h) += p25*q;
	  A(x+h,y-h) += p25*q;
	  A(x-h,y+h) += p25*q;
	  A(x+h,y+h) += p25*q;
	}
      }
    }

    /* Outermost 4 corners left unchanged. */

  } else if (job == 2) {

    /*******************************
     **                           **
     **  APPLY INVERSE TRANSFORM  **
     **                           **
     *******************************/

    /* Walk from smallest scales to largest ones. */
    for (s=2 ; s<=n ; s+=s) {

      /* Half-step size. */
      h = s/2;

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	A(x,0) -= (A(x-h,0) + A(x+h,0) + A(x,h))*p33;
      }
      for (y=h ; y<n ; y+=h) {
	A(0,y) -= (A(0,y-h) + A(0,y+h) + A(h,y))*p33;
	for (x=s ; x<n ; x+=s) {
	  A(x,y) -= (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
	}
	A(n,y) -= (A(n,y-h) + A(n,y+h) + A(n-h,y))*p33;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  A(x,y) -= (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
	}
      }
      for (x=h ; x<n ; x+=s) {
	A(x,n) -= (A(x,n-h) + A(x-h,n) + A(x+h,n))*p33;
      }

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  A(x,y) -= (A(x-h,y-h) + A(x+h,y-h) + A(x-h,y+h) + A(x+h,y+h))*p25;
	}
      }
    }

    /* Outermost 4 corners left unchanged. */

  } else if (job == 3) {

    /*****************************************
     **                                     **
     **  APPLY INVERSE TRANSPOSE TRANSFORM  **
     **                                     **
     *****************************************/

    /* Walk from largest scales to smallest ones. */
    for (s=n ; s>=2 ; s/=2) {

      /* Half-step size. */
      h = s/2;

      /* center of squares */
      for (y=h ; y<n ; y+=s) {
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y)*p25;
	  A(x-h,y-h) -= q;
	  A(x+h,y-h) -= q;
	  A(x-h,y+h) -= q;
	  A(x+h,y+h) -= q;
	}
      }

      /* centers of lozenges and borders */
      for (x=h ; x<n ; x+=s) {
	q = A(x,0)*p33;
	A(x-h,0) -= q;
	A(x+h,0) -= q;
	A( x ,h) -= q;
      }
      for (y=h ; y<n ; y+=h) {
	q = A(0,y)*p33;
	A(0,y-h) -= q;
	A(0,y+h) -= q;
	A(h,y)   -= q;
	for (x=s ; x<n ; x+=s) {
	  q = A(x,y)*p25;
	  A( x ,y-h) -= q;
	  A(x-h, y ) -= q;
	  A(x+h, y ) -= q;
	  A( x ,y+h) -= q;
	}
	q = A(n,y)*p33;
	A(n,y-h) -= q;
	A(n,y+h) -= q;
	A(n-h,y) -= q;
	if ((y += h) >= n) break;
	for (x=h ; x<n ; x+=s) {
	  q = A(x,y)*p25;
	  A( x ,y-h) -= q;
	  A(x-h, y ) -= q;
	  A(x+h, y ) -= q;
	  A( x ,y+h) -= q;
	}
      }
      for (x=h ; x<n ; x+=s) {
	q = A(x,n)*p33;
	A(x,n-h) -= q;
	A(x-h,n) -= q;
	A(x+h,n) -= q;
      }
    }

  } else {

    RETURN_FAILURE(BAD_JOB);

  }

#undef A

  RETURN_SUCCESS;
}
#endif /* FRIM_INT_2D */

/*
 * Return an array with structure function:
 *
 *    ALPHA*(R/R0)^BETA
 *
 * forward transform:
 *   DST is a DIM*DIM array
 *   SRC is a DIM*DIM + 2 array (two first coefficients are the normalized
 *                               tilts)
 *   DIM is the side dimension of the phase screen, must be a power
 *   of 2 plus 1.
 *
 *
 *   DIM = 2^P + 1
 *   NP = 2*P + 1        **** FIXME ****
 *
 *   D - array of structure funtion values for different separations.
 *   ND - Number of elements of D.
 *
 *   D[0] = SF(STP)
 *   D[1] = SF(sqrt(2)*STP)
 *   D[2] = SF(2*STP)
 *   ...
 *   D[k] = SF(STP*2^(k/2))
 *   ...
 *   D[ND - 2] = SF((DIM - 1)*STP)
 *   D[ND - 1] = SF(sqrt(2)*(DIM - 1)*STP)
 *
 * where STP is the pixel size, SF(r) is the structure function for a
 * separation r:
 *
 *   SF(r) = < [f(r' + r) - f(r')]^2 >
 *
 * for instance for a Kolmogorov phase screen:
 *
 *   SF(r) = 6.88 * (r/r0)^(5/3)
 *
 */

#ifdef FRIM_LANE
int FRIM_LANE(real_t *dst, const real_t *src, const size_t dim,
	      const real_t d[], size_t nd, unsigned int flags,
	      int *status)
{
  const real_t O_5  = 0.5;
  const real_t O_25 = 0.25;

  real_t d1, d2, d3, w1, w2, w3, p1, p2, p3, p4;
  size_t x, y, n, s, h, k, p;

#define SRC(i1,i2) src[(i2)*dim + (i1)]
#define DST(i1,i2) dst[(i2)*dim + (i1)]

  /* Check the arguments. */
  n = dim - 1;
  for (p = 0 ; 1<<p != n ; ) {
    if (++p >= (8*sizeof(size_t))) {
      RETURN_FAILURE(BAD_DIM);
    }
  }
  if (nd != 2*p + 2) {
    RETURN_FAILURE(BAD_NCOEFS);
  }

  /* Generate the 4 first phases. */
  k = nd;
  d2 = d[--k]; /* structure function at sqrt(2)*(DIM - 1)*STP */
  d1 = d[--k]; /* structure function at         (DIM - 1)*STP */
  if (flags & FRIM_SIX_VALUES) {
    /* Use 6 random values to generate the 4 first corners. */
    w1 = sqrt(d1 - d2/2.0);   /* standard deviation of corners */
    w2 = sqrt((d2 - d1)/2.0); /* standard deviation of diagonals */

    /* Get the tilts. */
    p1 = w1*(*src++);
    p2 = w1*(*src++);

    /* Generate the 4 first phases. */
    DST(0,0) = SRC(0,0)*w2 + p1;
    DST(n,n) = SRC(n,n)*w2 - p1;
    DST(n,0) = SRC(n,0)*w2 + p2;
    DST(0,n) = SRC(0,n)*w2 - p2;
  } else {
    /* Only use 4 random values. */
    real_t var = d1/2.0; /* variance such that corners on the same side
			    are uncorrelated */
    w1 = sqrt(var - d1/4.0 - d2/8.0);
    w2 = sqrt(d1/4.0 - d2/8.0);
    w3 = sqrt(d2/4.0);
    p1 = w1*SRC(0,0); /* piston */
    p2 = w2*SRC(n,0);
    p3 = w3*SRC(0,n); /* tip-tilt 1 */
    p4 = w3*SRC(n,n); /* tip-tilt 2 */
    DST(0,0) = p1 - p2 - p3;
    DST(n,0) = p1 + p2 - p4;
    DST(0,n) = p1 + p2 + p4;
    DST(n,n) = p1 - p2 + p3;
  }

  /* Generate other phases. */
  for (s=n ; s>=2 ; s/=2) {
    /* shift scale */
    h  = s/2;
    d3 = d2;
    d2 = d1;
    d1 = d[--k];
    w1 = sqrt(d1 - d2/4.0 - d3/8.0);

    /* centers of squares */
    for (y=h; y<n; y+=s) {
      for (x=h; x<n; x+=s) {
	DST(x,y) = SRC(x,y)*w1
	  + (DST(x-h, y-h) + DST(x-h, y+h) +
	     DST(x+h, y-h) + DST(x+h, y+h))*O_25;
      }
    }

    /* shift scale */
    d3 = d2;
    d2 = d1;
    d1 = d[--k];
    w1 = sqrt(d1 - d2/4.0 - d3/8.0);
    w2 = sqrt(d1 - d2/4.0);

    /* centers of lozenges */
    for (y=h; y < n; y+=s) {
      for (x=s; x<n; x+=s) {
	DST(x,y) = SRC(x,y)*w1
	  +(DST(x,y-h) + DST(x,y+h) + DST(x-h,y) + DST(x+h,y))*O_25;
      }
    }
    for (y=s; y < n; y+=s) {
      for (x=h; x<n; x+=s) {
	DST(x,y) = SRC(x,y)*w1
	  + (DST(x,y-h) +DST(x,y+h) + DST(x-h,y) +DST(x+h,y))*O_25;
      }
    }

    /* Borders. */
    for (y=h; y < n; y+=s) {
      DST(0,y) = SRC(0,y)*w2 + (DST(0,y-h) + DST(0,y+h))*O_5;
      DST(n,y) = SRC(n,y)*w2 + (DST(n,y-h) + DST(n,y+h))*O_5;
    }
    for (x=h; x<n; x+=s) {
      DST(x,0) = SRC(x,0)*w2 + (DST(x-h,0) + DST(x+h,0))*O_5;
      DST(x,n) = SRC(x,n)*w2 + (DST(x-h,n) + DST(x+h,n))*O_5;
    }
  }
  RETURN_SUCCESS;
}
#undef SRC
#undef DST
#endif /* FRIM_LANE */


#undef suffix
#undef real_t
#undef FRIM_INT_2D
#undef FRIM_GEN_2D
#undef FRIM_LANE

#endif /* _FRIM_CODE ----------------------------------------------------*/
