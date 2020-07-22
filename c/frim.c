/*
 * frim.c -
 *
 * FRiM (FRactal Iterative Method) main code.
 *
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of FRiM which is licensed under the MIT "Expat" License.
 *
 * Copyright (c) 2005-2020: Éric Thiébaut <https://github.com/emmt>
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
		if (status) *status = FRIM_PASTE(FRIM_ERROR_, CODE); \
		return FRIM_FAILURE
#define RETURN_SUCCESS \
		if (status) *status = FRIM_ERROR_NONE; \
		return FRIM_SUCCESS

#define TYPE    float
#define SUFFIX  _f
#define SQRT    sqrtf
#include __FILE__

#define TYPE    double
#define SUFFIX  _d
#define SQRT    sqrt
#include __FILE__

#else /* _FRIM_CODE defined ---------------------------------------------*/

#define WFS_MODEL1   FRIM_PASTE(wfs_model1, SUFFIX)
#define FRIM_GEN_2D  FRIM_PASTE(frim_gen_2d, SUFFIX)
#define FRIM_INT_2D  FRIM_PASTE(frim_int_2d, SUFFIX)
#define FRIM_LANE    FRIM_PASTE(frim_lane, SUFFIX)

#ifdef WFS_MODEL1
int WFS_MODEL1(TYPE*restrict dst, TYPE const*restrict src,
               size_t const dim, unsigned int flags, int* status)
{
    const TYPE p5 = 0.5;

    if (dim < 1) {
        RETURN_FAILURE(BAD_DIM);
    }
    size_t n = dim - 1;

    if ((flags & FRIM_JOB_TRANSPOSE)) {
#define SRC(s,i1,i2) src[((i2)*n + (i1))*2 + (s)]
#define DST(i1,i2)   dst[(i2)*dim + (i1)]
        if ((flags & FRIM_CLEAR)) {
            memset(dst, 0, dim*dim*sizeof(*dst));
        }
        for (size_t y = 0; y < n; ++y) {
            for (size_t x = 0; x < n; ++x) {
                TYPE wx = SRC(0,x,y)*p5;
                TYPE wy = SRC(1,x,y)*p5;
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
        for (size_t y = 0; y < n; ++y) {
            for (size_t x = 0; x < n; ++x) {
                DST(0,x,y) = (SRC(x+1,y+1) + SRC(x+1,y) - SRC(x,y+1) - SRC(x,y))*p5;
                DST(1,x,y) = (SRC(x+1,y+1) - SRC(x+1,y) + SRC(x,y+1) - SRC(x,y))*p5;
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
 *   C[2] = Cov(STP*√2)
 *   ...
 *   C[k] = Cov(STP*2^((k-1)/2))  for k = 1,..., NC - 1
 *   ...
 *   C[NC - 2] = Cov(STP*(DIM - 1))
 *   C[NC - 1] = Cov(STP*√2*(DIM - 1))
 *
 * NC - Number of elements for array of precomputed covariances.
 *
 * DIM = 2^P + 1
 * NC = 2*P + 3
 */

#ifdef FRIM_GEN_2D
int FRIM_GEN_2D(TYPE*restrict dst, TYPE const*restrict src,
                size_t const dim, TYPE const*restrict c,
                size_t const nc, int job, int* status)
{
    /* Check arguments. */
    if (dim < 1) {
        RETURN_FAILURE(BAD_DIM);
    }
    size_t n = dim - 1;
    size_t p = 0;
    while ((1 << p) < n && p < 8*sizeof(size_t) - 1) {
        ++p;
    }
    if ((1 << p) != n) {
        RETURN_FAILURE(BAD_DIM);
    }
    if (nc != 2*p + 3) {
        RETURN_FAILURE(BAD_NCOEFS);
    }

    /* Copy SRC in DST array to perform 'in-place' operation. */
    if (src && src != dst) {
        memcpy(dst, src, dim*dim*sizeof(src[0]));
    }

#define A(i1,i2) dst[(i2)*dim + (i1)]

    TYPE c0 = c[0]; /* C0 <- Cov(0) = Var */
    TYPE c1, c2, c3;
    TYPE w1, w2, w3, w4, w5;
    TYPE q;
    TYPE p1, p2, p3, p4;

    if (job == 0) {

        /*******************************
         **                           **
         **  APPLY FORWARD TRANSFORM  **
         **                           **
         *******************************/

        /* Walk from largest scales to smallest ones. */
        size_t k = nc;

        /* Outermost 4 corners. */
        c2 = c[--k]; /* C2 <- Cov(√2*(DIM - 1)) */
        c1 = c[--k]; /* C1 <- Cov(DIM - 1) */
        w1 = SQRT(c0 + 2*c1 + c2)/2;
        w2 = SQRT(c0 - 2*c1 + c2)/2;
        w3 = SQRT((c0 - c2)/2);
        p1 = w1*A(0,0);
        p2 = w2*A(n,0);
        p3 = w3*A(0,n);
        p4 = w3*A(n,n);
        A(0,0) = p1 - p2 - p3;
        A(n,0) = p1 + p2 - p4;
        A(0,n) = p1 + p2 + p4;
        A(n,n) = p1 - p2 + p3;

        for (size_t s = n; s >= 2; s /= 2) {
            /* Shift scale factors. */
            size_t h = s/2;
            c3 = c2;
            c2 = c1;
            c1 = c[--k];
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
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
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            q = c1/((c0 + c3)*c0 - 2*c2*c2);
            w3 = (c0 - c2)*q;
            w4 = (c0 - 2*c2 + c3)*q;
            w5 = SQRT(c0 - (w3 + w4 + w3)*c1);

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                A(x,0) = (A(x-h,0) + A(x+h,0))*w3 + A(x,h)*w4 + A(x,0)*w5;
            }
            for (size_t y = h; y < n; y += h) {
                A(0,y) = (A(0,y-h) + A(0,y+h))*w3 + A(h,y)*w4 + A(0,y)*w5;
                for (size_t x = s; x < n; x += s) {
                    A(x,y) = (A( x ,y-h) +
                              A(x-h, y ) +
                              A(x+h, y ) +
                              A( x ,y+h))*w1 + A(x,y)*w2;
                }
                A(n,y) = (A(n,y-h) + A(n,y+h))*w3 + A(n-h,y)*w4 + A(n,y)*w5;
                if ((y += h) >= n) break;
                for (size_t x = h; x < n; x += s) {
                    A(x,y) = (A( x ,y-h) +
                              A(x-h, y ) +
                              A(x+h, y ) +
                              A( x ,y+h))*w1 + A(x,y)*w2;
                }
            }
            for (size_t x = h; x < n; x += s) {
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
        size_t k = 0;
        c2 = c[++k];
        c3 = c[++k];

        for (size_t s = 2; s <= n; s += s) {

            /* Half-step size. */
            size_t h = s/2;

            /* Shift scale factors. */
            c1 = c2;
            c2 = c3;
            c3 = c[++k];
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            q = c1/((c0 + c3)*c0 - 2*c2*c2);
            w3 = (c0 - c2)*q;
            w4 = (c0 - 2*c2 + c3)*q;
            w5 = SQRT(c0 - (w3 + w4 + w3)*c1);

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                q = A(x,0);
                A( x ,0)  = w5*q;
                A(x-h,0) += w3*q;
                A(x+h,0) += w3*q;
                A( x ,h) += w4*q;
            }
            for (size_t y = h; y < n; y += h) {
                q = A(0,y);
                A(0,y)    = w5*q;
                A(0,y-h) += w3*q;
                A(0,y+h) += w3*q;
                A(h,y)   += w4*q;
                for (size_t x = s; x < n; x += s) {
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
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y);
                    A( x,  y )  = w2*q;
                    A( x ,y-h) += w1*q;
                    A(x-h, y ) += w1*q;
                    A(x+h, y ) += w1*q;
                    A( x ,y+h) += w1*q;
                }
            }
            for (size_t x = h; x < n; x += s) {
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
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
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
        c2 = c3; /* C2 <- Cov(√2*(DIM - 1)) */
        w1 = SQRT(c0 + 2*c1 + c2)*0.5;
        w2 = SQRT(c0 - 2*c1 + c2)*0.5;
        w3 = SQRT((c0 - c2)*0.5);
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
        size_t k = 0;
        c2 = c[++k];
        c3 = c[++k];

        for (size_t s = 2; s <= n; s += s) {

            /* Half-step size. */
            size_t h = s/2;

            /* Shift scale factors. */
            c1 = c2;
            c2 = c3;
            c3 = c[++k];
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            w1 /= w2;
            w2 = 1/w2;
            q = c1/((c0 + c3)*c0 - 2*c2*c2);
            w3 = (c0 - c2)*q;
            w4 = (c0 - 2*c2 + c3)*q;
            w5 = SQRT(c0 - (w3 + w4 + w3)*c1);
            w3 /= w5;
            w4 /= w5;
            w5 = 1/w5;

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                A(x,0) = A(x,0)*w5 - (A(x-h,0) + A(x+h,0))*w3 - A(x,h)*w4;
            }
            for (size_t y = h; y < n; y += h) {
                A(0,y) = A(0,y)*w5 - (A(0,y-h) + A(0,y+h))*w3 - A(h,y)*w4;
                for (size_t x = s; x < n; x += s) {
                    A(x,y) = A(x,y)*w2 - (A( x ,y-h) +
                                          A(x-h, y ) +
                                          A(x+h, y ) +
                                          A( x ,y+h))*w1;
                }
                A(n,y) = A(n,y)*w5 - (A(n,y-h) + A(n,y+h))*w3 - A(n-h,y)*w4;
                if ((y += h) >= n) break;
                for (size_t x = h; x < n; x += s) {
                    A(x,y) = A(x,y)*w2 - (A( x ,y-h) +
                                          A(x-h, y ) +
                                          A(x+h, y ) +
                                          A( x ,y+h))*w1;
                }
            }
            for (size_t x = h; x < n; x += s) {
                A(x,n) = A(x,n)*w5 - A(x,n-h)*w4 - (A(x-h,n) + A(x+h,n))*w3;
            }

            /* Shift scale factors. */
            c1 = c2;
            c2 = c3;
            c3 = c[++k];
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            w1 /= w2;
            w2 = 1/w2;

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
                    A(x,y) = A(x,y)*w2 - (A(x-h,y-h) +
                                          A(x+h,y-h) +
                                          A(x-h,y+h) +
                                          A(x+h,y+h))*w1;
                }
            }
        }

        /* Outermost 4 corners. */
        c1 = c2; /* C1 <- Cov(DIM - 1) */
        c2 = c3; /* C2 <- Cov(√2*(DIM - 1)) */
#if 0 /* unused because seems to produce larger rounding errors */
        w1 = SQRT(0.25/(c0 + 2*c1 + c2));
        w2 = SQRT(0.25/(c0 - 2*c1 + c2));
        w3 = SQRT(0.5/(c0 - c2));
#else
        w1 = 1/SQRT((c0 + 2*c1 + c2)*4);
        w2 = 1/SQRT((c0 - 2*c1 + c2)*4);
        w3 = 1/SQRT((c0 - c2)*2);
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
        size_t k = nc;

        /* Outermost 4 corners. */
        c2 = c[--k]; /* C2 <- Cov(√2*(DIM - 1)) */
        c1 = c[--k]; /* C1 <- Cov(DIM - 1) */
#if 0
        w1 = SQRT((c0 + 2*c1 + c2)*4);
        w2 = SQRT((c0 - 2*c1 + c2)*4);
#endif
        p1 = A(0,0)/SQRT((c0 + 2*c1 + c2)*4);
        p2 = A(n,0)/SQRT((c0 - 2*c1 + c2)*4);
        w3 = SQRT((c0 - c2)*2);
        p3 = A(0,n)/w3;
        p4 = A(n,n)/w3;
        A(0,0) = p1 - p2 - p3;
        A(n,0) = p1 + p2 - p4;
        A(0,n) = p1 + p2 + p4;
        A(n,n) = p1 - p2 + p3;

        for (size_t s = n; s >= 2; s /= 2) {
            /* Shift scale factors. */
            size_t h = s/2;
            c3 = c2;
            c2 = c1;
            c1 = c[--k];
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            w1 /= w2;
            w2 = 1/w2;

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
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
            w1 = c1/(c0 + 2*c2 + c3);
            w2 = SQRT(c0 - 4*c1*w1);
            w1 /= w2;
            w2 = 1/w2;
            q = c1/((c0 + c3)*c0 - 2*c2*c2);
            w3 = (c0 - c2)*q;
            w4 = (c0 - 2*c2 + c3)*q;
            w5 = SQRT(c0 - (w3 + w4 + w3)*c1);
            w3 /= w5;
            w4 /= w5;
            w5 = 1/w5;

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                q = A(x,0);
                A( x ,0)  = w5*q;
                A(x-h,0) -= w3*q;
                A(x+h,0) -= w3*q;
                A( x ,h) -= w4*q;
            }
            for (size_t y = h; y < n; y += h) {
                q = A(0,y);
                A(0,y) = w5*q;
                A(0,y-h) -= w3*q;
                A(0,y+h) -= w3*q;
                A(h,y)   -= w4*q;
                for (size_t x = s; x < n; x += s) {
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
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y);
                    A(x,y) = w2*q;
                    A( x ,y-h) -= w1*q;
                    A(x-h, y ) -= w1*q;
                    A(x+h, y ) -= w1*q;
                    A( x ,y+h) -= w1*q;
                }
            }
            for (size_t x = h; x < n; x += s) {
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
int FRIM_INT_2D(TYPE*restrict dst, TYPE const*restrict src,
                size_t const dim, int job, int* status)
{

    /* Check arguments. */
    if (dim < 1) {
        RETURN_FAILURE(BAD_DIM);
    }
    size_t n = dim - 1;
    size_t p = 0;
    while ((1 << p) < n && p < 8*sizeof(size_t) - 1) {
        ++p;
    }
    if ((1 << p) != n) {
        RETURN_FAILURE(BAD_DIM);
    }

    /* Copy SRC in DST array to perform 'in-place' operation. */
    if (src && src != dst) {
        memcpy(dst, src, dim*dim*sizeof(src[0]));
    }

#define A(i1,i2) dst[(i2)*dim + (i1)]

    const TYPE p25 = 1.0/4.0;
    const TYPE p33 = 1.0/3.0;
    TYPE q;

    if (job == 0) {

        /*******************************
         **                           **
         **  APPLY FORWARD TRANSFORM  **
         **                           **
         *******************************/

        /* Walk from largest scales to smallest ones.
           Outermost 4 corners are left unchanged. */
        for (size_t s = n; s >= 2; s /= 2) {

            /* Half-step size. */
            size_t h = s/2;

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
                    A(x,y) += (A(x-h,y-h) + A(x+h,y-h) + A(x-h,y+h) + A(x+h,y+h))*p25;
                }
            }

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                A(x,0) += (A(x-h,0) + A(x+h,0) + A(x,h))*p33;
            }
            for (size_t y = h; y < n; y += h) {
                A(0,y) += (A(0,y-h) + A(0,y+h) + A(h,y))*p33;
                for (size_t x = s; x < n; x += s) {
                    A(x,y) += (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
                }
                A(n,y) += (A(n,y-h) + A(n,y+h) + A(n-h,y))*p33;
                if ((y += h) >= n) break;
                for (size_t x = h; x < n; x += s) {
                    A(x,y) += (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
                }
            }
            for (size_t x = h; x < n; x += s) {
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
        for (size_t s = 2; s <= n; s += s) {

            /* Half-step size. */
            size_t h = s/2;

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                q = A(x,0)*p33;
                A(x-h,0) += q;
                A(x+h,0) += q;
                A( x ,h) += q;
            }
            for (size_t y = h; y < n; y += h) {
                q = A(0,y)*p33;
                A(0,y-h) += q;
                A(0,y+h) += q;
                A(h,y)   += q;
                for (size_t x = s; x < n; x += s) {
                    q = A(x,y)*p25;
                    A( x ,y-h) += q;
                    A(x-h, y ) += q;
                    A(x+h, y ) += q;
                    A( x ,y+h) += q;
                }
                q = A(n,y)*p33;
                A( n ,y-h) += q;
                A( n ,y+h) += q;
                A(n-h, y ) += q;
                if ((y += h) >= n) break;
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y)*p25;
                    A( x ,y-h) += q;
                    A(x-h, y ) += q;
                    A(x+h, y ) += q;
                    A( x ,y+h) += q;
                }
            }
            for (size_t x = h; x < n; x += s) {
                q = A(x,n)*p33;
                A( x ,n-h) += q;
                A(x-h, n ) += q;
                A(x+h, n ) += q;
            }

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y)*p25;
                    A(x-h,y-h) += q;
                    A(x+h,y-h) += q;
                    A(x-h,y+h) += q;
                    A(x+h,y+h) += q;
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
        for (size_t s = 2; s <= n; s += s) {

            /* Half-step size. */
            size_t h = s/2;

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                A(x,0) -= (A(x-h,0) + A(x+h,0) + A(x,h))*p33;
            }
            for (size_t y = h; y < n; y += h) {
                A(0,y) -= (A(0,y-h) + A(0,y+h) + A(h,y))*p33;
                for (size_t x = s; x < n; x += s) {
                    A(x,y) -= (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
                }
                A(n,y) -= (A(n,y-h) + A(n,y+h) + A(n-h,y))*p33;
                if ((y += h) >= n) break;
                for (size_t x = h; x < n; x += s) {
                    A(x,y) -= (A(x,y-h) + A(x-h,y) + A(x+h,y) + A(x,y+h))*p25;
                }
            }
            for (size_t x = h; x < n; x += s) {
                A(x,n) -= (A(x,n-h) + A(x-h,n) + A(x+h,n))*p33;
            }

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
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
        for (size_t s = n; s >= 2; s /= 2) {

            /* Half-step size. */
            size_t h = s/2;

            /* center of squares */
            for (size_t y = h; y < n; y += s) {
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y)*p25;
                    A(x-h,y-h) -= q;
                    A(x+h,y-h) -= q;
                    A(x-h,y+h) -= q;
                    A(x+h,y+h) -= q;
                }
            }

            /* centers of lozenges and borders */
            for (size_t x = h; x < n; x += s) {
                q = A(x,0)*p33;
                A(x-h,0) -= q;
                A(x+h,0) -= q;
                A( x ,h) -= q;
            }
            for (size_t y = h; y < n; y += h) {
                q = A(0,y)*p33;
                A(0,y-h) -= q;
                A(0,y+h) -= q;
                A(h,y)   -= q;
                for (size_t x = s; x < n; x += s) {
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
                for (size_t x = h; x < n; x += s) {
                    q = A(x,y)*p25;
                    A( x ,y-h) -= q;
                    A(x-h, y ) -= q;
                    A(x+h, y ) -= q;
                    A( x ,y+h) -= q;
                }
            }
            for (size_t x = h; x < n; x += s) {
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
 *   D[1] = SF(√2*STP)
 *   D[2] = SF(2*STP)
 *   ...
 *   D[k] = SF(STP*2^(k/2))
 *   ...
 *   D[ND - 2] = SF((DIM - 1)*STP)
 *   D[ND - 1] = SF(√2*(DIM - 1)*STP)
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
int FRIM_LANE(TYPE*restrict dst, TYPE const*restrict src,
              size_t const dim, TYPE const*restrict d, size_t nd,
              unsigned int flags, int *status)
{
    /* Check arguments. */
    if (dim < 1) {
        RETURN_FAILURE(BAD_DIM);
    }
    size_t n = dim - 1;
    size_t p = 0;
    while ((1 << p) < n && p < 8*sizeof(size_t) - 1) {
        ++p;
    }
    if ((1 << p) != n) {
        RETURN_FAILURE(BAD_DIM);
    }
    if (nd != 2*p + 2) {
        RETURN_FAILURE(BAD_NCOEFS);
    }

#define SRC(i1,i2) src[(i2)*dim + (i1)]
#define DST(i1,i2) dst[(i2)*dim + (i1)]

    const TYPE p5  = 0.5;
    const TYPE p25 = 0.25;
    TYPE d1, d2, d3, w1, w2, w3, p1, p2, p3, p4;

    /* Generate the 4 first phases. */
    size_t k = nd;
    d2 = d[--k]; /* structure function at √2*(DIM - 1)*STP */
    d1 = d[--k]; /* structure function at    (DIM - 1)*STP */
    if (flags & FRIM_SIX_VALUES) {
        /* Use 6 random values to generate the 4 first corners. */
        w1 = SQRT(d1 - d2/2);   /* standard deviation of corners */
        w2 = SQRT((d2 - d1)/2); /* standard deviation of diagonals */

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
        TYPE var = d1/2; /* variance such that corners on the same side
                              are uncorrelated */
        w1 = SQRT(var - d1/4 - d2/8);
        w2 = SQRT(d1/4 - d2/8);
        w3 = SQRT(d2/4);
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
    for (size_t s = n; s >= 2; s /= 2) {
        /* shift scale */
        size_t h  = s/2;
        d3 = d2;
        d2 = d1;
        d1 = d[--k];
        w1 = SQRT(d1 - d2/4 - d3/8);

        /* centers of squares */
        for (size_t y = h; y < n; y += s) {
            for (size_t x = h; x < n; x += s) {
                DST(x,y) = SRC(x,y)*w1
                    + (DST(x-h, y-h) + DST(x-h, y+h) +
                       DST(x+h, y-h) + DST(x+h, y+h))*p25;
            }
        }

        /* shift scale */
        d3 = d2;
        d2 = d1;
        d1 = d[--k];
        w1 = SQRT(d1 - d2/4 - d3/8);
        w2 = SQRT(d1 - d2/4);

        /* centers of lozenges */
        for (size_t y = h; y < n; y += s) {
            for (size_t x = s; x < n; x += s) {
                DST(x,y) = SRC(x,y)*w1
                    +(DST(x,y-h) + DST(x,y+h) + DST(x-h,y) + DST(x+h,y))*p25;
            }
        }
        for (size_t y = s; y < n; y += s) {
            for (size_t x = h; x < n; x += s) {
                DST(x,y) = SRC(x,y)*w1
                    + (DST(x,y-h) +DST(x,y+h) + DST(x-h,y) +DST(x+h,y))*p25;
            }
        }

        /* Borders. */
        for (size_t y = h; y < n; y += s) {
            DST(0,y) = SRC(0,y)*w2 + (DST(0,y-h) + DST(0,y+h))*p5;
            DST(n,y) = SRC(n,y)*w2 + (DST(n,y-h) + DST(n,y+h))*p5;
        }
        for (size_t x = h; x < n; x += s) {
            DST(x,0) = SRC(x,0)*w2 + (DST(x-h,0) + DST(x+h,0))*p5;
            DST(x,n) = SRC(x,n)*w2 + (DST(x-h,n) + DST(x+h,n))*p5;
        }
    }

#undef SRC
#undef DST

    RETURN_SUCCESS;
}
#endif /* FRIM_LANE */


#undef SUFFIX
#undef TYPE
#undef SQRT
#undef FRIM_INT_2D
#undef FRIM_GEN_2D
#undef FRIM_LANE

#endif /* _FRIM_CODE ----------------------------------------------------*/
