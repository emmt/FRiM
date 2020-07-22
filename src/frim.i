/*
 * frim.i -
 *
 * FRiM (FRactal Iterative Method) in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of FRiM which is licensed under the MIT "Expat" License.
 *
 * Copyright (c) 2005-2020: Éric Thiébaut <https://github.com/emmt>
 *
 *-----------------------------------------------------------------------------
 */

if (is_func(plug_in)) plug_in, "frim";

/*---------------------------------------------------------------------------*/

if (! is_func(random_n)) {
  include, "random.i";
}

func kolmogorov_lane(p, r0=, six=)
/* DOCUMENT kolmogorov_lane(p);

     Generate a (2^P + 1)-by-(2^P + 1) Kolmogorov phase screen by Lane et
     al. (1992) mid-point algorithm.  Keyword R0 can be used to set the
     Fried's parameter in pixels (default: R0=1).  Keyword SIX has the same
     meaning as in frim_lane (which see).

   SEE ALSO: frim_lane, structure_function.
 */
{
    /* Generate structure function for a Kolmogorov phase screen. */
    dim = 2^p + 1;
    alpha = 6.88;
    beta = 5.0/3.0;
    if (! is_void(r0)) alpha *= (1.0/r0)^beta;
    return frim_lane(alpha*sqrt(2.0)^(beta*indgen(0:2*p + 1)), six=six);
}

func frim_lane(sf, six=)
/* DOCUMENT frim_lane(sf, six=);

     Generate a 2-D array of random values obeying the structure function
     SF by Lane's et al. [1] mid-point algorithm.  The result is a
     DIM-by-DIM array where DIM = 2^P + 1 and where P is such that
     numberof(SF) = 2*(P + 1).  The structure function SF of a function F
     is defined as:

          SF(r) = < (f(r' + r) - f(r'))^2 >

     where < > denotes averaging over r'.  SF must be given for:

          R = [1,sqrt(2),2,..., sqrt(2)*(DIM-1)]
            = sqrt(2.0)^indgen(0:2*p+1)

     If keyword SIX is true, original method is used which requires 6
     random values to generate the 4 values at the outermost 4 corner of
     the result; otherwise, a modified method is used which only requires
     as many random values as the result array.


   REFERENCES:
     [1]  R. G. Lane, A. lindemann and J. C. Dainty, "Simulation of a
          Kolmogorov phase screen", in "Wave in random media", vol. 2,
          pp. 209--224, 1992.

   SEE ALSO: random_n, structure_function.
 */
{
    if ((type = structof(sf)) == double) {
        op = _frim_lane_d;
    } else if (type == float || type==long || type==int ||
               type==short || type==char) {
        type = float;
        op = _frim_lane_f;
    } else {
        error, "bad data type for structure function";
    }
    n = numberof(sf);
    p = n/2 - 1;
    if (p < 1 || 2*p + 2 != n) {
        error, "bad number of elements for structure function";
    }
    dim = (1 << p) + 1;
    res = array(double, dim, dim);
    if (six) {
        rnd = random_n(2 + dim*dim);
        flags = 2n;
    } else {
        rnd = random_n(dim, dim);
        flags = 0n;
    }
    status = 0n;
    if (op(res, rnd, dim, sf, n, flags, status)) {
        frim_error, status;
    }
    return res;
}

extern _frim_lane_f;
/* PROTOTYPE
   int frim_lane_f(float array dst, float array src, long dim,
                    float array d, long nd, int flags,
                    int array status);
*/
extern _frim_lane_d;
/* PROTOTYPE
   int frim_lane_d(double array dst, double array src, long dim,
                    double array d, long nd, int flags,
                    int array status);
*/

/*---------------------------------------------------------------------------*/

func kolmogorov_covariance(p, r0=)
/* DOCUMENT kolmogorov_covariance(p);

     Returns the spatial covariance COV(r) of a Kolmogorov phase screen PHI
     for distances:

       R = R0*sqrt(2.0)^indgen(0 : 2*p + 2)

     Keyword R0 can be used to set the Fried's parameter in pixels
     (default: R0=1).  The covariance is defined as:

       COV(r) = < PHI(r' + r)*PHI(r') >
              = COV(0) - 1/2 SF(r)

     where < > denotes averaging over r' and where SF(r) is the structure
     function of the phase PHI.  The returned array can be used to generate
     a random phase screen with a uniform variance (unlike Lane's et
     al. method) thanks to the function frim_gen_2d:

        frim_gen_2d(COV, random_n(DIM, DIM))

     with DIM = 2^P + 1.


   SEE ALSO: kolmogorov_lane, frim_gen_2d, random_n.
 */
{
    alpha = 6.88;
    beta = 5.0/3.0;
    if (! is_void(r0)) alpha *= (1.0/r0)^beta;

    /* Generate covariances for a Kolmogorov phase screen. */

    /* Compute: W(r) = 3.44*(r/r0)^(5/3)
     * so that the covariances write: C(r) = C(0) - W(r)
     * W(r) is 1/2 times the structure function
     */
    nc = 2*p + 3;
    w = array(double, nc);
    if (nc > 1) {
        q = alpha/2;
        w(2) = q;
        if (nc > 2) w(3:nc) = q*2.0^((beta/2)*indgen(nc - 2));
    }

    /* The 'piston' variance is:
     *
     *   C_piston = C0 + 2*C2 + C3
     *            = 4*C0 - 2*W((DIM - 1)) - W(sqrt(2)*(DIM - 1))
     *
     * we could choose C0 so that the variance of the piston-like mode
     * is zero (i.e. no 'piston'):
     *
     *   C0 = 0.5*W((DIM - 1)) + 0.25*W(sqrt(2)*(DIM - 1))
     *
     * but this yields a non-invertible covariance matrix.  We therefore
     * define C0 so that all covariances are positive:
     */
    c0 = max(w);

    return c0 - w;
}

func kolmogorov_random(p, r0=)
/* DOCUMENT kolmogorov_random(p);

     Generate a (2^P + 1)-by-(2^P + 1) Kolmogorov phase screen with uniform
     variance by a fractal method.  Keyword R0 can be used to set the
     Fried's parameter in pixels (default: R0=1).

   SEE ALSO: frim_gen_2d_random, kolmogorov_covariance,
             frim_lane, structure_function.
 */
{
    return frim_gen_2d_random(kolmogorov_covariance(p, r0=r0));
}

func frim_gen_2d_random(cov)
/* DOCUMENT frim_gen_2d_random(cov);

     Generate a random 2-D field following the given covariance COV which
     must have 2*P + 3 elements.  The result is a DIM-by-DIM array with
     DIM = 2^P + 1.  See frim_gen_2d() for a description of the argument COV.

   SEE ALSO: frim_int_2d, kolmogorov_covariance, random_n.
 */
{
    dim = (1 << ((numberof(cov) - 3)/2)) + 1;
    if ((type = structof(cov)) == double) {
        op = _frim_gen_2d_d;
        x = random_n(dim, dim);
    } else if (type == float || type == long || type == int || type == short ||
               type == char) {
        op = _frim_gen_2d_f;
        x = float(random_n(dim, dim));
    } else {
        error, "bad data type for COV";
    }
    status = int();
    if (op(x, x, dim, cov, numberof(cov), 0, status)) frim_error, status;
    return x;
}

func frim_new_2d_linop(cov)
/* DOCUMENT frim_new_2d_linop();
         or frim_new_2d_linop(cov);

    Returns new linear operator (see doc. for linop_new) wrapped around
    frim_int_2d, if COV is unspecified; or around frim_gen_2d with
    covariance given by COV, otherwise.

   SEE ALSO: linop_new, frim_int_2d, frim_gen_2d.
 */
{
    if (is_void(cov)) {
        return linop_new(frim_int_2d);
    } else {
        return linop_new(frim_gen_2d, cov);
    }
}

func frim_gen_2d(cov, x, job)
/* DOCUMENT frim_gen_2d(cov, x);
         or frim_gen_2d(cov, x, job);

     Apply fractal interpolation operator to 2-D array X with given
     covariance COV.  X must be a DIM-by-DIM array with DIM = 2^P + 1 and
     COV must have 2*P + 3 elements.  The covariance of a field F is
     defined as:

         COV(r) = < F(r' + r)*F(r') >
                = COV(0) - 1/2 SF(r)

     where < > denotes averaging over r' and where SF(r) is the structure
     function of the field F.  The input covariance must be given for:

         R = sqrt(2)^indgen(0 : 2*P + 2)

     Argument JOB defines which variant of the linear operator is applied:

         JOB=0 or unspecified  - direct fractal interpolation
             1                 - transpose operator
             2                 - inverse operator
             3                 - inverse transpose operator

     For instance, with JOB=0 (or unspecified) and with X a DIM-by-DIM,
     array of pseudo-random values with normal distribution (see random_n),
     the result is a random 2-D field with covariance (approximatively)
     equals to COV.  See kolmogorov_covariance as an example for computing
     COV.

     Another usage of this function is to compute the regularization term
     in inverse problems for which the a priori covariance of the unknown
     is COV:

         REGUL(x) = sum(frim_gen_2d(COV, x, 2)^2);

      of which the gradient w.r.t. X is:

         GRAD(x) = 2*frim_gen_2d(COV, frim_gen_2d(COV, x, 2), 3)


   SEE ALSO: frim_int_2d, frim_gen_2d_random, kolmogorov_covariance, random_n.
 */
{
    if (is_void(job)) job = 0n;
    if (! is_array(x) || (dims = dimsof(x))(1) != 2 ||
        dims(2) != dims(3)) error, "expecting a N×N array for X";
    dim = dims(3);
    if ((type = structof(x)) == double) {
        op = _frim_gen_2d_d;
        x = double(x); /* make a private copy for in-place operation */
    } else if (type == float || type == long || type == int || type == short ||
               type == char) {
        op = _frim_gen_2d_f;
        x = float(x); /* make a private copy for in-place operation */
    } else {
        error, "bad data type for X";
    }
    status = int();
    if (op(x, x, dim, cov, numberof(cov), job, status)) frim_error, status;
    return x;
}

extern _frim_gen_2d_f;
/* PROTOTYPE
   int frim_gen_2d_f(float array dst, float array src,
                     long dim, float array c, long nc,
                     int job, int array status);
*/
extern _frim_gen_2d_d;
/* PROTOTYPE
   int frim_gen_2d_d(double array dst, double array src,
                     long dim, double array c, long nc,
                     int job, int array status);
*/

/*---------------------------------------------------------------------------*/

func frim_int_2d(x, job)
/* DOCUMENT frim_int_2d(x);
         or frim_int_2d(x, job);

     Apply fractal interpolation operator to 2-D array X.  X must be a
     DIM-by-DIM array with DIM = 2^P + 1.  Argument JOB defines which
     variant of the linear operator is applied:

         JOB=0 or unspecifed   - direct fractal interpolation
             1                 - transpose operator
             2                 - inverse operator
             3                 - inverse transpose operator


   SEE ALSO: frim_gen_2d.
 */
{
    if (is_void(job)) job = 0n;
    if (! is_array(x) || (dims = dimsof(x))(1) != 2 ||
        dims(2) != dims(3)) error, "expecting a N×N array for X";
    dim = dims(3);
    if ((type = structof(x)) == double) {
        op = _frim_int_2d_d;
        x = double(x); /* make a private copy for in-place operation */
    } else if (type == float || type == long || type == int || type == short ||
               type == char) {
        op = _frim_int_2d_f;
        x = float(x); /* make a private copy for in-place operation */
    } else {
        error, "bad data type for X";
    }
    status = int();
    if (op(x, x, dim, job, status)) frim_error, status;
    return x;
}

extern _frim_int_2d_f;
/* PROTOTYPE
   int frim_int_2d_f(float array dst, float array src,
                     long dim, int job, int array status);
*/
extern _frim_int_2d_d;
/* PROTOTYPE
   int frim_int_2d_d(double array dst, double array src,
                     long dim, int job, int array status);
*/

/*---------------------------------------------------------------------------*/

func frim_region_2d(dim)
/* DOCUMENT frim_region_2d(dim);
     Return a hash table (require Yeti) as follows:
        reg = frim_region_2d(dim);
        reg.np is an integer DIM-by-DIM array such that reg.np(x,y) is the
               number of parent cells used by FRIM to generate field at (x,y)
        reg.rp is an real DIM-by-DIM array such that reg.rp(x,y) is the
               distance to the parent cells

   SEE ALSO: h_new, frim_gen_2d.
 */
{
    p = long(log(dim -1)/log(2) + 0.5);
    if (dim != (1<<p) + 1) error, "bad dimension";
    n = dim - 1;
    np = array(long, dim, dim);
    rp = array(double, dim, dim);

    /* S = step size, or current scale
     * H = S/2 = distance of neighbors
     */
    for (s = n; s >= 2; s /= 2) {
        write, s;
        h = s/2;
        u = indgen(h+1:n:s);

        /* center of squares */
        np(u, u) = 4;
        rp(u, u) = sqrt(2)*h;

        /* borders */
        np(u, 1) = 3;
        rp(u, 1) = h;
        np(u, dim) = 3;
        rp(u, dim) = h;
        np(1, u) = 3;
        rp(1, u) = h;
        np(dim, u) = 3;
        rp(dim, u) = h;

        /* center of lozenges */
        if (s < n) {
            v = indgen(1+s:n:s);
            np(u, v) = 4;
            rp(u, v) = h;
            np(v, u) = 4;
            rp(v, u) = h;
        }
    }

    /* N = number of parent cells
     * R = distance to parents
     */
    return h_new(n=np, r=rp);
}

/* ---------------------------------------------------------------------- */
/* STRUCTURE FUNCTION */

/* Structure function of function F for a step S:
 *
 *   SF(s) = NUM(s)/DEN(s)
 *
 * where:
 *
 *   NUM(s) = sum_{r} w(r) w(r+s) [f(r) - f(r+s)]^2
 *   DEN(s) = sum_{r} w(r) w(r+s)
 *
 * W(R) is the support of F(R), DEN(S) is therefore a normalization term
 * equals to the number of pairs in the support for a separation S.
 *
 * Using Fourier transforms (IFT is the inverse Fourier transform):
 *   DEN(s) = IFT(conj(FT(w))*FT(w))
 *          = IFT(|FT(w)|^2)
 *   NUM(s) = sum_{r} w(r) w(r+s) [f(r) - f(r+s)]^2
 *          = sum_{r} [w(r+s) w(r) f(r)^2 + w(r) w(r+s) f(r+s)^2
 *                     - 2 w(r) f(r) w(r+s) f(r+s)]
 *          = IFT(FT(w)*conj(FT(w*f^2)) + conj(FT(w))*FT(w*f^2)
 *                - 2*FT(w*f)*conj(FT(w*f)))
 *          = 2*IFT(Re(FT(w)*conj(FT(w*f^2))) - |FT(w*f)|^2)
 *
 * because (intercorrelation and convolution theorems):
 *   sum_r conj(f(r)) g(r + s) = IFT(conj(FT(f))*FT(g))
 *   sum_r f(r) g(r - s) = IFT(FT(f)*FT(g))
 *
 */
func structure_function(f, support=, extra=, number=, integrate=)
/* DOCUMENT structure_function(f);

     Compute structure function of field F.  F can be a 1-D or 2-D real
     valued array or a function (or the name of a function) which returns a
     random field realization.  Keyword SUPPORT can be used to specify a
     support shape (SUPPORT must be conformable with fields and should be
     made of 1's and 0's).  Note that the result is rolled.

     If F is a function, keyword NUMBER can be set to the number of random
     fields to generate and the average structure function is returned.
     Keyword EXTRA can be used to specify any argument for F, i.e. a random
     field is generated by: F(EXTRA).

     If keyword INTEGRATE is true, then the array [SUM_USF, SUM_NRM] is
     returned where SUM_USF is the sum of the 'unnormalized' structure
     functions computed for all the generated fields and SUN_NRM is the
     normalization term (proportional to NUMBER).  The 'normalized'
     structure function can be obtained by: SUM_USF/SUM_NRM.  For instance,
     if two calls to wfs_strf are made for some number of random fields
     (not necessarily the sames):
        T1 = structure_function(F, extra=ARG, number=N1, integrate=1);
        T2 = structure_function(F, extra=ARG, number=N2, integrate=1);
     then the average structure function is:
        (T1(..,1) + T2(..,1))/(T1(..,2) + T2(..,2))

   SEE ALSO: fftw.
 */
{
    /* Initialization and default settings. */
    local wf;
    if (is_void(number)) number = 1;

    /* Get the first field. */
    if (structof(f) == string) f = symbol_def(f);
    if (is_func(f)) {
        fn = f;
        f = fn(extra);
    } else if (number != 1) {
        error, "F must be a function with NUMBER != 1";
    }
    if (! is_array(f)) {
        error, "non-array F";
    }
    dims = dimsof(f);
    ndims = dims(1);
    if (ndims == 1) {
        n1 = dims(2);
        l1 = fft_best_dim(2*n1 - 1);
        rng1 = 1:n1;
        rng2 = [];
        dimlist = [1,l1];
    } else if (ndims == 2) {
        n1 = dims(2);
        n2 = dims(3);
        l1 = fft_best_dim(2*n1 - 1);
        l2 = fft_best_dim(2*n2 - 1);
        rng1 = 1:n1;
        rng2 = 1:n2;
        dimlist = [2,l1,l2];
    } else {
        error, "too many dimensions";
    }

    /* Workspace and FFTW plan for forward real->complex FFT */
    tmp = array(double, dimlist);
    plan = fftw_plan(dimlist, +1, real=1);

    /* Loop over random generation. */
    nevals = 1;
    tmp(rng1, rng2) = (is_void(support) ? 1.0 : support);
    z0 = fftw(tmp, plan); /* FFT of W (only needed once) */
    z0_re = double(z0);
    z0_im = z0.im;
    z0 = []; /* free memory */
    for (;;) {
        if (is_void(support)) eq_nocopy, wf, f; else wf = support*f;
        tmp(rng1, rng2) = wf;
        z1 = fftw(tmp, plan); /* FFT of W*F */
        z1_re = double(z1);
        z1_im = z1.im;
        z1 = []; /* free memory */
        tmp(rng1, rng2) = wf*f;
        z2 = fftw(tmp, plan); /* FFT of W*F^2 */
        z2_re = double(z2);
        z2_im = z2.im;
        z2 = []; /* free memory */
        if (nevals == 1) {
            num = z0_re*z2_re + z0_im*z2_im - z1_re*z1_re - z1_im*z1_im;
        } else {
            num += (z0_re*z2_re + z0_im*z2_im - z1_re*z1_re - z1_im*z1_im);
        }

        if (nevals >= number) break;
        f = fn(extra);
        ++nevals;
    }

    /* FFTW plan for backward complex->real FFT */
    plan = fftw_plan(dimlist, -1, real=1);

    num = 2.0*fftw(num, plan);
    den = number*fftw(z0_re*z0_re + z0_im*z0_im, plan);

    /* Extract relevant part. */
    if (ndims >= 1) {
        if (l1 > 2*n1 - 1) {
            i = where(abs(fft_indgen(l1)) < n1);
            num = num(i,..);
            den = den(i,..);
        }
        if (ndims >= 2) {
            if (l2 > 2*n2 - 1) {
                i = where(abs(fft_indgen(l2)) < n2);
                num = num(,i,..);
                den = den(,i,..);
            }
        }
    }

    /* Return result (avoiding division by zero). */
    if (integrate) return [num, den];
    amax = (is_void(support) ? 1.0 : max(max(support), -min(support)));
    threshold = 1e-4*number*amax^2;
    if (min(den) >= threshold) return num/den;
    i = where(den >= threshold);
    (scl = array(double, dimsof(den)))(i) = 1.0/den(i);
    return scl*num;
}

/*---------------------------------------------------------------------------*/
/* STREHL RATIO */

func strehl_ratio(pherr, pupil=, zap_chessboard=)
/* DOCUMENT strehl_ratio(pherr);

     Return Strehl ratio for phase error PHERR, computed as:

         exp(-avg((PHERR - avg(PHERR))^2))

     Keyword PUPIL can be set with an array of same geometry of PHERR and
     which is true (non-zero) inside the pupil and false elsewhere.

     If keyword ZAP_CHESSBOARD is true, the 'chessboard' modes are removed
     from PHERR.

   SEE ALSO: kolmogorov_lane, phase_variance.
 */
{
    /* account for pupil */
    if (! is_void(pupil)) pherr = pherr(where(pupil));

    /* remove piston */
    pherr -= avg(pherr);
    if (zap_chessboard) {
        /* remove chessboard modes */
        w = array(0.5, dimsof(pherr));
        w(where(indgen(numberof(w)) % 2)) = -0.5;
        pherr -= (sum(w*pherr)/sum(w*w))*w;
        w = 1.0 - w;
        pherr -= (sum(w*pherr)/sum(w*w))*w;
    }
    return exp(-avg(pherr*pherr));
}

func phase_variance(pherr, pupil=, zap_chessboard=)
/* DOCUMENT phase_variance(pherr);

     Return phase variance for phase error PHERR, computed as:

         avg((PHERR - avg(PHERR))^2)

     Keyword PUPIL can be set with an array of same geometry of PHERR and
     which is true (non-zero) inside the pupil and false elsewhere.

     If keyword ZAP_CHESSBOARD is true, the 'chessboard' modes are removed
     from PHERR.

   SEE ALSO: kolmogorov_lane, strehl_ratio.
 */
{
    /* account for pupil */
    if (! is_void(pupil)) pherr = pherr(where(pupil));

    /* remove piston */
    pherr -= avg(pherr);
    if (zap_chessboard) {
        /* remove chessboard modes */
        w = array(0.5, dimsof(pherr));
        w(where(indgen(numberof(w)) % 2)) = -0.5;
        pherr -= (sum(w*pherr)/sum(w*w))*w;
        w = 1.0 - w;
        pherr -= (sum(w*pherr)/sum(w*w))*w;
    }
    return avg(pherr*pherr);
}

/*---------------------------------------------------------------------------*/
/* WAVEFRONT SENSOR MODEL */

func wfs_model1(x, tr)
/* DOCUMENT wfs_model1(p);
         or wfs_model1(r, 1);

     Apply wavefront sensor model (No. 1) operator. P is a N-by-N array and
     the result is a 2-by-(N-1)-by-(N-1) array Y such that:
        Y(1,..) = X(dif,zcen)
        Y(2,..) = X(zcen,dif)
     If the second argument is true, the transpose operator is applied, R
     must be a 2-by-(N-1)-by-(N-1) array and the result is a N-by-N array.

   SEE ALSO:
 */
{
    if ((type = structof(x)) == double) {
        op = _wfs_model1_d;
    } else {
        op = _wfs_model1_f;
        if (type != float) {
            if (type==char || type==short || type==int || type==long) {
                type = float;
                x = float(x);
            } else {
                error, "bad data type for X";
            }
        }
    }
    if (tr) {
        if (! is_array(x) || (dims = dimsof(x))(1) != 3 || dims(2) != 2 ||
            dims(3) != dims(4)) error, "expecting a 2×N×N array for X";
        dim = dims(3) + 1;
        flags = 1n;
        y = array(type, dim, dim);
    } else {
        if (! is_array(x) || (dims = dimsof(x))(1) != 2 ||
            dims(2) != dims(3)) error, "expecting a N×N array for X";
        dim = dims(3);
        flags = 0n;
        y = array(type, 2, dim - 1, dim - 1);
    }
    status = int();
    if (op(y, x, dim, flags, status)) frim_error, status;
    return y;
}

extern _wfs_model1_f;
/* PROTOTYPE
   int wfs_model1_f(float array dst, float array src, long dim,
                    int flags, int array status);
*/
extern _wfs_model1_d;
/* PROTOTYPE
   int wfs_model1_d(double array dst, double array src, long dim,
                    int flags, int array status);
*/

/*---------------------------------------------------------------------------*/
extern frim_error;
/* DOCUMENT frim_error, code;

     Generate error from FRIM package with message according to CODE.

   SEE ALSO: error.
*/
