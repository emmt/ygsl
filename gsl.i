/*
 * gsl.i --
 *
 * Support for GSL (GNU Scientific Library) in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2012 Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license as
 * circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 */

if (is_func(plug_in)) plug_in, "ygsl";

local gsl_sf;
/* DOCUMENT gsl_sf_*
 *
 *   Special functions from GSL (GNU Scientific Library) are prefixed with
 *   "gsl_sf_"; to obtain more information, see the following documentation
 *   entries:
 *
 *     gsl_sf_airy_Ai   - Airy functions
 *     gsl_sf_bessel_J0 - regular cylindrical Bessel functions
 *     gsl_sf_bessel_Y0 - irregular cylindrical Bessel functions
 *     gsl_sf_bessel_I0 - regular modified cylindrical Bessel functions
 *     gsl_sf_bessel_K0 - irregular modified cylindrical Bessel functions
 *     gsl_sf_bessel_j0 - regular spherical Bessel functions
 *     gsl_sf_bessel_y0 - irregular spherical Bessel functions
 *     gsl_sf_bessel_i0_scaled - regular modified spherical Bessel functions
 *     gsl_sf_bessel_k0_scaled - irregular modified spherical Bessel functions
 *     gsl_sf_clausen - Clausen function
 *     gsl_sf_dawson - Dawson integral
 *     gsl_sf_debye - Debye functions
 *     gsl_sf_dilog - dilogarithm
 *     gsl_sf_ellint_Kcomp - Legendre form of complete elliptic integrals
 *     gsl_sf_erf - error functions
 *     gsl_sf_exp - exponential and logarithm functions
 *     gsl_sf_expint - exponential, hyperbolic and trigonometric integrals
 *     gsl_sf_fermi_dirac - Fermi-Dirac integrals
 *     gsl_sf_gamma - Gamma functions
 *     gsl_sf_psi - Digamma, trigamma and polygamma functions
 *     gsl_sf_lamber - Lambert's functions
 *     gsl_sf_legendre - Legendre polynomials
 *     gsl_sf_synchrotron - synchrotron functions
 *     gsl_sf_transport - transport functions
 *     gsl_sf_sin - trigonometric functions
 *     gsl_sf_zeta - Zeta functions
 */

extern gsl_sf_airy_Ai;
extern gsl_sf_airy_Bi;
extern gsl_sf_airy_Ai_scaled;
extern gsl_sf_airy_Bi_scaled;
extern gsl_sf_airy_Ai_deriv;
extern gsl_sf_airy_Bi_deriv;
extern gsl_sf_airy_Ai_deriv_scaled;
extern gsl_sf_airy_Bi_deriv_scaled;
/* DOCUMENT gsl_sf_airy_Ai(x [,flags])
 *          gsl_sf_airy_Bi(x [,flags])
 *          gsl_sf_airy_Ai_deriv(x [,flags])
 *          gsl_sf_airy_Bi_deriv(x [,flags])
 *          gsl_sf_airy_Ai_scaled(x [,flags])
 *          gsl_sf_airy_Bi_scaled(x [,flags])
 *          gsl_sf_airy_Ai_deriv_scaled(x [,flags])
 *          gsl_sf_airy_Bi_deriv_scaled(x [,flags])
 *
 *   These routines compute the Airy functions and derivatives for the
 *   argument X (a non-complex numerical array).
 *
 *   The routines gsl_sf_airy_Ai and gsl_sf_airy_Bi compute Airy functions
 *   Ai(x) and Bi(x) which are defined by the integral representations:
 *
 *      Ai(x) = (1/PI) \int_0^\infty cos((1/3)*t^3 + x*t) dt
 *      Bi(x) = (1/PI) \int_0^\infty (exp(-(1/3)*t^3)
 *                                    + sin((1/3)*t^3 + x*t)) dt
 *
 *   The routines gsl_sf_airy_Ai_deriv and gsl_sf_airy_Bi_deriv compute
 *   the derivatives of the Airy functions.
 *
 *   The routines gsl_sf_airy_Ai_scaled and gsl_sf_airy_Bi_scaled compute
 *   a scaled version of the Airy functions S_A(x) Ai(x) and S_B(x) Bi(x).
 *   The scaling factors are:
 *      S_A(x) = exp(+(2/3)*x^(3/2)), for x>0
 *               1,                   for x<0;
 *      S_B(x) = exp(-(2/3)*x^(3/2)), for x>0
 *               1,                   for x<0.
 *
 *   The routines gsl_sf_airy_Ai_deriv_scaled and
 *   gsl_sf_airy_Bi_deriv_scaled compute the derivatives of the scaled Airy
 *   functions.
 *
 *   The optional FLAGS argument is a bitwise combination which specifies
 *   the relative accuracy of the result and if an estimate of the error
 *   is required:
 *
 *     (FLAGS & 1) is non-zero to compute an estimate of the error, the
 *         result, says Y, has an additional dimension of length 2
 *         prepended to the dimension list of X:
 *             Y(1,..) = value of F(X)
 *             Y(2,..) = error estimate for the value of F(X)
 *
 *     (FLAGS & 6) is the accuracy mode:
 *         6 - Double-precision (GSL_PREC_DOUBLE), a relative accuracy of
 *             approximately 2e-16.
 *         4 - Single-precision (GSL_PREC_SINGLE), a relative accuracy of
 *             approximately 1e-7.
 *         2 - Approximate values (GSL_PREC_APPROX), a relative accuracy
 *             of approximately 5e-4.
 *         0 - Default accuracy (GSL_PREC_DOUBLE).
 *
 *   For instance, with FLAGS=1, function values are computed with relative
 *   accuracy of 2e-16 and an estimate of the error is returned; with
 *   FLAGS=2, approximate values with relative accuracy of 5e-4 are
 *   returned without error estimate
 *
 *
 * SEE ALSO: gsl_sf.
 */

extern gsl_sf_bessel_J0;
extern gsl_sf_bessel_J1;
extern gsl_sf_bessel_Jn;
extern gsl_sf_bessel_Jnu;
/* DOCUMENT gsl_sf_bessel_J0(x [,err])
 *          gsl_sf_bessel_J1(x [,err])
 *          gsl_sf_bessel_Jn(n, x [,err])
 *          gsl_sf_bessel_Jnu(nu, x [,err])
 *
 *   These functions compute the regular cylindrical Bessel functions for
 *   argument X (a non-complex numerical array or scalar) and of various
 *   order: zeroth order, J_0(x); first order, J_1(x), integer order order
 *   N, J_n(x), and fractional order NU, J_nu(x).  N must be a scalar
 *   integer and NU a scalar real.
 *
 *   If optional argument ERR is true, these functions also compute an
 *   estimate of the error, the result, says Y, has an additional dimension
 *   of length 2 prepended to the dimension list of X:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_Y0, gsl_sf_bessel_I0, gsl_sf_bessel_K0,
 *           gsl_sf_bessel_j0, gsl_sf_bessel_y0, gsl_sf_bessel_i0,
 *           gsl_sf_bessel_k0.
 */

extern gsl_sf_bessel_Y0;
extern gsl_sf_bessel_Y1;
extern gsl_sf_bessel_Yn;
extern gsl_sf_bessel_Ynu;
/* DOCUMENT gsl_sf_bessel_Y0(x [,err])
 *          gsl_sf_bessel_Y1(x [,err])
 *          gsl_sf_bessel_Yn(n, x [,err])
 *          gsl_sf_bessel_Ynu(nu, x [,err])
 *
 *   These functions compute the irregular cylindrical Bessel functions for
 *   X>0.  See gsl_sf_bessel_J0 for a more detailled description of the
 *   arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_I0;
extern gsl_sf_bessel_I1;
extern gsl_sf_bessel_In;
extern gsl_sf_bessel_Inu;
extern gsl_sf_bessel_I0_scaled;
extern gsl_sf_bessel_I1_scaled;
extern gsl_sf_bessel_In_scaled;
extern gsl_sf_bessel_Inu_scaled;
/* DOCUMENT gsl_sf_bessel_I0(x [,err])
 *          gsl_sf_bessel_I1(x [,err])
 *          gsl_sf_bessel_In(n, x [,err])
 *          gsl_sf_bessel_Inu(nu, x [,err])
 *          gsl_sf_bessel_I0_scaled(x [,err])
 *          gsl_sf_bessel_I1_scaled(x [,err])
 *          gsl_sf_bessel_In_scaled(n, x [,err])
 *          gsl_sf_bessel_Inu_scaled(nu, x [,err])
 *
 *   These routines compute the regular modified cylindrical Bessel
 *   functions and their scaled counterparts.  The scaling factor is
 *   exp(-abs(X)); for instance: I0_scaled(X) = exp(-abs(X))*I0(X).  See
 *   gsl_sf_bessel_J0 for a more detailled description of the arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_K0;
extern gsl_sf_bessel_K1;
extern gsl_sf_bessel_Kn;
extern gsl_sf_bessel_Knu;
extern gsl_sf_bessel_lnKnu;
extern gsl_sf_bessel_K0_scaled;
extern gsl_sf_bessel_K1_scaled;
extern gsl_sf_bessel_Kn_scaled;
extern gsl_sf_bessel_Knu_scaled;
/* DOCUMENT gsl_sf_bessel_K0(x [,err])
 *          gsl_sf_bessel_K1(x [,err])
 *          gsl_sf_bessel_Kn(n, x [,err])
 *          gsl_sf_bessel_Knu(nu, x [,err])
 *          gsl_sf_bessel_lnKnu(nu, x [,err])
 *          gsl_sf_bessel_K0_scaled(x [,err])
 *          gsl_sf_bessel_K1_scaled(x [,err])
 *          gsl_sf_bessel_Kn_scaled(n, x [,err])
 *          gsl_sf_bessel_Knu_scaled(nu, x [,err])
 *
 *   These routines compute the irregular modified cylindrical Bessel
 *   functions and their scaled counterparts.  The scaling factor is exp(X)
 *   for X>0; for instance: K0_scaled(X) = exp(X)*K0(X).  The function
 *   gsl_sf_bessel_lnKnu computes the logarithm of the irregular modified
 *   Bessel function of fractional order NU.  See gsl_sf_bessel_J0 for a
 *   more detailled description of the arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_j0;
extern gsl_sf_bessel_j1;
extern gsl_sf_bessel_j2;
extern gsl_sf_bessel_jl;
/* DOCUMENT gsl_sf_bessel_j0(x [,err])
 *          gsl_sf_bessel_j1(x [,err])
 *          gsl_sf_bessel_j2(x [,err])
 *          gsl_sf_bessel_jl(l, x [,err])
 *
 *   These routines compute the regular spherical Bessel functions of
 *   zeroth order (j0), first order (j1), second order (j2) and l-th order
 *   (jl, for X>=0 and L>=0).  See gsl_sf_bessel_J0 for a more detailled
 *   description of the arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_y0;
extern gsl_sf_bessel_y1;
extern gsl_sf_bessel_y2;
extern gsl_sf_bessel_yl;
/* DOCUMENT gsl_sf_bessel_y0(x [,err])
 *          gsl_sf_bessel_y1(x [,err])
 *          gsl_sf_bessel_y2(x [,err])
 *          gsl_sf_bessel_yl(l, x [,err])
 *
 *   These routines compute the irregular spherical Bessel functions of
 *   zeroth order (y0), first order (y1), second order (y2) and l-th order
 *   (yl, for L>=0):
 *
 *     y0(x) = -cos(x)/x
 *     y1(x) = -(cos(x)/x + sin(x))/x
 *     y2(x) = (-3/x^3 + 1/x)*cos(x) - (3/x^2)*sin(x)
 *
 *   See gsl_sf_bessel_J0 for a more detailled description of the
 *   arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_i0_scaled;
extern gsl_sf_bessel_i1_scaled;
extern gsl_sf_bessel_i2_scaled;
extern gsl_sf_bessel_il_scaled;
/* DOCUMENT gsl_sf_bessel_i0_scaled(x [,err])
 *          gsl_sf_bessel_i1_scaled(x [,err])
 *          gsl_sf_bessel_i2_scaled(x [,err])
 *          gsl_sf_bessel_il_scaled(l, x [,err])
 *
 *   These routines compute the regular modified spherical Bessel functions
 *   of zeroth order (i0), first order (i1), second order (i2) and l-th
 *   order (il):
 *
 *     il_scaled(x) = exp(-abs(x))*il(x)
 *
 *   The regular modified spherical Bessel functions i_l(x) are related to
 *   the modified Bessel functions of fractional order by:
 *
 *     i_l(x) = sqrt(PI/(2*x))*I_{l + 1/2}(x)
 *
 *   See gsl_sf_bessel_J0 for a more detailled description of the
 *   arguments.
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_bessel_k0_scaled;
extern gsl_sf_bessel_k1_scaled;
extern gsl_sf_bessel_k2_scaled;
extern gsl_sf_bessel_kl_scaled;
/* DOCUMENT gsl_sf_bessel_k0_scaled(x [,err])
 *          gsl_sf_bessel_k1_scaled(x [,err])
 *          gsl_sf_bessel_k2_scaled(x [,err])
 *          gsl_sf_bessel_kl_scaled(l, x [,err])
 *
 *   These routines compute the irregular modified spherical Bessel
 *   functions of zeroth order (k0), first order (k1), second order (k2)
 *   and l-th order (kl), for X>0:
 *
 *     kl_scaled(x) = exp(x)*kl(x)
 *
 *   The irregular modified spherical Bessel functions i_l(x) are related to
 *   the modified Bessel functions of fractional order by:
 *
 *     k_l(x) = sqrt(PI/(2*x))*K_{l + 1/2}(x)
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_bessel_J0.
 */

extern gsl_sf_clausen;
/* DOCUMENT gsl_sf_clausen(x [,err])
 *
 *   Returns the Clausen function Cl_2 of its argument X.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

extern gsl_sf_dawson;
/* DOCUMENT gsl_sf_dawson(x [,err])
 *
 *   Returns the Dawson integral of its argument X defined by:
 *
 *       exp(-x^2) \int_0^x exp(t^2) dt
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

extern gsl_sf_debye_1;
extern gsl_sf_debye_2;
extern gsl_sf_debye_3;
extern gsl_sf_debye_4;
extern gsl_sf_debye_5;
extern gsl_sf_debye_6;
local gsl_sf_debye;
/* DOCUMENT gsl_sf_debye_1(x [,err])
 *          gsl_sf_debye_2(x [,err])
 *          gsl_sf_debye_3(x [,err])
 *          gsl_sf_debye_4(x [,err])
 *          gsl_sf_debye_5(x [,err])
 *          gsl_sf_debye_6(x [,err])
 *
 *   Return the Debye function D_n(x) of argument X defined by the
 *   following integral:
 *
 *     D_n(x) = n/x^n \int_0^x (t^n/(e^t - 1)) dt
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

extern gsl_sf_dilog;
/* DOCUMENT gsl_sf_dilog(x [,err])
 *
 *   Return the dilogarithm for a real argument X.  If optional argument
 *   ERR is true, the result, says Y, has an additional dimension of length
 *   2 prepended to the dimension list of X which is used to provide an
 *   estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

extern gsl_sf_ellint_Kcomp;
extern gsl_sf_ellint_Ecomp;
/* DOCUMENT gsl_sf_ellint_Kcomp(k [,flags])
 *          gsl_sf_ellint_Ecomp(k [,flags])
 *   Return the complete elliptic integral K(k) or E(k).  See
 *   gsl_sf_airy_Ai for the meaning of optional argument FLAGS.
 *
 * SEE ALSO: gsl_sf, gsl_sf_airy_Ai.
 */

extern gsl_sf_erf;
extern gsl_sf_erfc;
extern gsl_sf_log_erfc;
extern gsl_sf_erf_Z;
extern gsl_sf_erf_Q;
extern gsl_sf_hazard;
/* DOCUMENT gsl_sf_erf(x [,err])
 *          gsl_sf_erfc(x [,err])
 *          gsl_sf_log_erfc(x [,err])
 *          gsl_sf_erf_Q(x [,err])
 *          gsl_sf_erf_Z(x [,err])
 *          gsl_sf_hazard(x [,err])
 *
 *   gsl_sf_erf(x) computes the error function:
 *
 *       erf(x) = (2/sqrt(pi)) \int_0^x exp(-t^2) dt
 *
 *   gsl_sf_erfc(x) computes the complementary error function:
 *
 *       erfc(x) = 1 - erf(x)
 *               = (2/sqrt(pi)) \int_x^\infty exp(-t^2) dt
 *
 *   gsl_sf_log_erfc(x) computes the logarithm of the complementary error function.
 *
 *   gsl_sf_erf_Z(x) computes the Gaussian probability density function:
 *
 *       Z(x) = (1/sqrt(2 pi)) \exp(-x^2/2).
 *
 *   gsl_sf_erf_Q(x) computes the upper tail of the Gaussian probability
 *   density function:
 *
 *       Q(x) = (1/sqrt(2 pi)) \int_x^\infty \exp(-t^2/2) dt.
 *
 *   gsl_sf_hazard(x) computes the hazard function for the normal
 *   distribution, also known as the inverse Mill's ratio:
 *
 *       h(x) = Z(x)/Q(x)
 *            = sqrt(2/pi) exp(-x^2/2)/erfc(x/sqrt(2)).
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

extern gsl_sf_exp;
extern gsl_sf_expm1;
extern gsl_sf_exprel;
extern gsl_sf_exprel_2;
extern gsl_sf_exprel_n;
extern gsl_sf_log;
extern gsl_sf_log_abs;
extern gsl_sf_log_1plusx;
extern gsl_sf_log_1plusx_mx;
/* DOCUMENT gsl_sf_exp(x [,err])
 *          gsl_sf_expm1(x [,err])
 *          gsl_sf_exprel(x [,err])
 *          gsl_sf_exprel_2(x [,err])
 *          gsl_sf_exprel_n(n, x [,err])
 *          gsl_sf_log(x [,err])
 *          gsl_sf_log_abs(x [,err])
 *          gsl_sf_log_1plusx(x [,err])
 *          gsl_sf_log_1plusx_mx(x [,err])
 *
 *   gsl_sf_exp(X) computes the exponential of X.
 *
 *   gsl_sf_expm1(X) computes the quantity exp(X) - 1 using an algorithm
 *   that is accurate for small X.
 *
 *   gsl_sf_exprel(X) computes the quantity (exp(X) - 1)/X using an
 *   algorithm that is accurate for small X and which is based on the
 *   expansion:
 *
 *       (exp(x) - 1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 *   gsl_sf_exprel_2(X) computes the quantity 2*(exp(X) - 1)/X^2 using an
 *   algorithm that is accurate for small X and which is based on the
 *   expansion:
 *
 *       2*(exp(x) - 1 - x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 *   gsl_sf_exprel_n(N,X) computes the N-relative exponential (N must be a
 *   scalar integer):
 *
 *       expre_n(x) = n! / x^n ( exp(x) - \sum_{k=0}^{n-1} x^k / k! )
 *
 *   gsl_sf_log(X) computes the logarithm of X, for X > 0.
 *
 *   gsl_sf_log_abs(X) computes the logarithm of |X|, for X != 0.
 *
 *   gsl_sf_log_1plusx(x) computes log(1 + X) for X > -1 using an algorithm
 *   that is accurate for small X.
 *
 *   gsl_sf_log_1plusx_mx(x) computes log(1 + X) - X for X > -1 using an
 *   algorithm that is accurate for small X.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

local gsl_sf_expint;
extern gsl_sf_expint_E1;
extern gsl_sf_expint_E2;
extern gsl_sf_expint_Ei;
extern gsl_sf_expint_3;
extern gsl_sf_Shi;
extern gsl_sf_Chi;
extern gsl_sf_Si;
extern gsl_sf_Ci;
extern gsl_sf_atanint;
/* DOCUMENT gsl_sf_expint_E1(x [, err])
 *          gsl_sf_expint_E2(x [, err])
 *          gsl_sf_expint_Ei(x [, err])
 *          gsl_sf_expint_3(x [, err])
 *          gsl_sf_Shi(x [, err])
 *          gsl_sf_Chi(x [, err])
 *          gsl_sf_Si(x [, err])
 *          gsl_sf_Ci(x [, err])
 *          gsl_sf_atanint(x [, err])
 *
 *   gsl_sf_expint_E1(X) computes the exponential integral:
 *       E1(x) = \int_1^\infty exp(-x t)/t dt
 *
 *   gsl_sf_expint_E2(X) computes the second-order exponential integral:
 *       E2(x) = \int_1^\infty exp(-x t)/t^2 dt
 *
 *   gsl_sf_expint_E2(X) computes the exponetial integral:
 *       Ei(x) = -PV( \int_{-x}^\infty exp(-t)/t dt )
 *   where PV() denotes the principal value.
 *
 *   gsl_sf_expint_3(X) computes the third-order exponential integral:
 *       Ei_3(x) = \int_0^x \exp(-t^3) dt       for x >= 0.
 *
 *   gsl_sf_Shi(X) computes the integral:
 *       Shi(x) = \int_0^x sinh(t)/t dt.
 *
 *   gsl_sf_Chi(X) computes the integral:
 *       Chi(x) = Re[ gamma_E + log(x) + \int_0^x (cosh(t) - 1)/t dt ]
 *   where gamma_E is the Euler constant.
 *
 *   gsl_sf_Si(X) computes the Sine integral:
 *       Si(x) = \int_0^x sin(t)/t dt.
 *
 *   gsl_sf_Ci(X) computes the Cosine integral:
 *       Ci(x) = -\int_x^\int_x cos(t)/t dt        for x > 0.
 *
 *   gsl_sf_atanint(X) computes the arc-tangent integral:
 *       AtanInt(x) = \int_0^x arctan(t)/t dt.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf
 */

local gsl_sf_fermi_dirac;
extern gsl_sf_fermi_dirac_m1;
extern gsl_sf_fermi_dirac_0;
extern gsl_sf_fermi_dirac_1;
extern gsl_sf_fermi_dirac_2;
extern gsl_sf_fermi_dirac_mhalf;
extern gsl_sf_fermi_dirac_half;
extern gsl_sf_fermi_dirac_3half;
extern gsl_sf_fermi_dirac_int;
/* DOCUMENT gsl_sf_fermi_dirac_int(j, x [, err])
 *          gsl_sf_fermi_dirac_m1(x [, err])
 *          gsl_sf_fermi_dirac_0(x [, err])
 *          gsl_sf_fermi_dirac_1(x [, err])
 *          gsl_sf_fermi_dirac_2(x [, err])
 *          gsl_sf_fermi_dirac_mhalf(x [, err])
 *          gsl_sf_fermi_dirac_half(x [, err])
 *          gsl_sf_fermi_dirac_3half(x [, err])
 *
 *   gsl_sf_fermi_dirac_int(J,X) computes the complete Fermi-Dirac integral
 *   with an index of J:
 *       F_j(x) = 1/Gamma(j + 1) \int_0^\infty t^j/(exp(t - x) + 1) dt
 *   where J is a scalar integer and Gamma() is the Gamma function:
 *       Gamma(n) = (n - 1)!
 *   for integer n.
 *
 *   gsl_sf_fermi_dirac_m1(X) computes the complete Fermi-Dirac integral
 *   with an index of -1:
 *       F_{-1}(x) = exp(x)/(1 + exp(x))
 *
 *   gsl_sf_fermi_dirac_0(X) computes the complete Fermi-Dirac integral
 *   with an index of 0:
 *       F_0(x) = log(1 + exp(x))
 *
 *   gsl_sf_fermi_dirac_1(X) computes the complete Fermi-Dirac integral
 *   with an index of 1:
 *       F_1(x) = \int_0^\infty t/(exp(t - x) + 1) dt
 *
 *   gsl_sf_fermi_dirac_2(X) computes the complete Fermi-Dirac integral
 *   with an index of 2:
 *       F_2(x) = (1/2) \int_0^\infty t^2/(exp(t - x) + 1) dt
 *
 *   gsl_sf_fermi_dirac_mhalf(X) computes the complete Fermi-Dirac integral
 *   with an index of -1/2.
 *
 *   gsl_sf_fermi_dirac_half(X) computes the complete Fermi-Dirac integral
 *   with an index of +1/2.
 *
 *   gsl_sf_fermi_dirac_3half(X) computes the complete Fermi-Dirac integral
 *   with an index of +3/2.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf, gsl_sf_gamma.
 */

extern gsl_sf_gamma;
extern gsl_sf_lngamma;
extern gsl_sf_gammastar;
extern gsl_sf_gammainv;
extern gsl_sf_taylorcoeff;
/* DOCUMENT gsl_sf_gamma(x [, err])
 *          gsl_sf_lngamma(x [, err])
 *          gsl_sf_gammastar(x [, err])
 *          gsl_sf_gammainv(x [, err])
 *          gsl_sf_taylorcoeff(n, x [, err])
 *
 *   gsl_sf_gamma(X) computes the Gamma function:
 *       Gammma(x) = \int_0^\infty t^(x - 1) exp(-t) dt          for x >= 0
 *   for a positive integer argument, Gamma(n) = (n - 1)!.
 *
 *   gsl_sf_lngamma(X) computes the logarithm of the Gamma function.
 *
 *   gsl_sf_gammastar(X) computes the regulated Gamma function:
 *       GammaStar(x) = Gamma(x) / ( sqrt(2 pi) x^(x - 1/2) exp(x) )
 *                    = 1 + 1/12x + ...     for large x
 *
 *   gsl_sf_gammainv(X) computes the reciprocal of the Gamma function
 *   1/Gamma(x) using the real Lanczos method.
 *
 *   gsl_sf_taylorcoeff(N,X) computes the Taylor coefficient X^N/N!
 *   for X >= 0 and N >= 0 -- N must be a scalar integer.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

extern gsl_sf_psi;
extern gsl_sf_psi_1piy;
extern gsl_sf_psi_1;
extern gsl_sf_psi_n;
/* DOCUMENT gsl_sf_psi(x [, err])
 *          gsl_sf_psi_1(x [, err])
 *          gsl_sf_psi_n(n, x [, err])
 *          gsl_sf_psi_1piy(x [, err])
 *
 *   gsl_sf_psi(X) computes the digamma function \psi(x) for X != 0.
 *
 *   gsl_sf_psi_1piy(Y) computes the real part of the digamma function on the
 *   line 1+i y, \Re[\psi(1 + i y)].
 *
 *   gsl_sf_psi_1(X) computes the trigamma function \psi'(x) for X.
 *
 *   gsl_sf_psi_n(N, X) computes the polygamma function \psi^{(n)}(x)
 *   for N >= 0, X > 0.
 *
 * SEE ALSO: gsl_sf.
 */

local gsl_sf_lambert;
extern gsl_sf_lambert_W0;
extern gsl_sf_lambert_Wm1;
/* DOCUMENT gsl_sf_lambert_W0(x [, err])
 *          gsl_sf_lambert_Wm1(x [, err])
 *   Lambert's W functions, W(x), are defined to be solutions of the
 *   equation W(x) exp(W(x)) = x.  This function has multiple branches for
 *   x < 0; however, it has only two real-valued branches.  We define W0(x)
 *   to be the principal branch, where W > -1 for x < 0, and Wm1(x) to
 *   be the other real branch, where W < -1 for x < 0.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

local gsl_sf_legendre;
extern gsl_sf_legendre_P1;
extern gsl_sf_legendre_P2;
extern gsl_sf_legendre_P3;
extern gsl_sf_legendre_Pl;
extern gsl_sf_legendre_Q0;
extern gsl_sf_legendre_Q1;
extern gsl_sf_legendre_Ql;
/* DOCUMENT gsl_sf_legendre_P1(x [, err])
 *          gsl_sf_legendre_P2(x [, err])
 *          gsl_sf_legendre_P3(x [, err])
 *          gsl_sf_legendre_Pl(l, x [, err])
 *          gsl_sf_legendre_Q0(x [, err])
 *          gsl_sf_legendre_Q1(x [, err])
 *          gsl_sf_legendre_Ql(l, x [, err])
 *
 *   The functions gsl_sf_legendre_P# evaluate the Legendre polynomials
 *   P_l(x) for specific values of l = 1, 2, 3 or for a scalar integer l.
 *
 *   The functions gsl_sf_legendre_Q# evaluate the Legendre function
 *   Q_l(x) for specific values of l = 0, 1 or for a scalar integer l.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

local gsl_sf_synchrotron;
extern gsl_sf_synchrotron_1;
extern gsl_sf_synchrotron_2;
local gsl_sf_transport;
extern gsl_sf_transport_2;
extern gsl_sf_transport_3;
extern gsl_sf_transport_4;
extern gsl_sf_transport_5;
/* DOCUMENT gsl_sf_synchrotron_1(x [, err])
 *          gsl_sf_synchrotron_2(x [, err])
 *          gsl_sf_transport_2(x [, err])
 *          gsl_sf_transport_3(x [, err])
 *          gsl_sf_transport_4(x [, err])
 *          gsl_sf_transport_5(x [, err])
 *
 *   gsl_sf_synchrotron_1(x) computes the first synchrotron function:
 *       x \int_x^\infty K_{5/3}(t) dt        for x >= 0.
 *
 *   gsl_sf_synchrotron_2(x) computes the second synchrotron function:
 *       x K_{2/3}(x)                         for x >= 0.
 *
 *   The transport functions J(n,x) are defined by the integral representations:
 *       J(n,x) = \int_0^x t^n e^t /(e^t - 1)^2 dt.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

extern gsl_sf_sin;
extern gsl_sf_cos;
extern gsl_sf_sinc;
extern gsl_sf_lnsinh;
extern gsl_sf_lncosh;
/* DOCUMENT gsl_sf_sin(x [, err])
 *          gsl_sf_cos(x [, err])
 *          gsl_sf_sinc(x [, err])
 *          gsl_sf_lnsinh(x [, err])
 *          gsl_sf_lncosh(x [, err])
 *
 *   gsl_sf_sin(X) computes the sine function of X.
 *
 *   gsl_sf_cos(X) computes the cosine function of X.
 *
 *   gsl_sf_sinc(X) computes sinc(x) = sin(pi x)/(pi x) for any value of X.
 *
 *   gsl_sf_lnsinh(X) computes log(sinh(X)) for X > 0.
 *
 *   gsl_sf_lncosh(X) computes log(cosh(X)) for any value of X.
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

extern gsl_sf_zeta;
extern gsl_sf_zetam1;
extern gsl_sf_eta;
/* DOCUMENT gsl_sf_zeta(x [, err])
 *          gsl_sf_zetam1(x [, err])
 *          gsl_sf_eta(x [, err])
 *
 *   gsl_sf_zeta(x) computes the Riemann zeta function:
 *       zeta(x) = \sum_{k=1}^\infty k^{-x}    for X != 1.
 *
 *   gsl_sf_zetam1(x) computes zeta(X) - 1 for X != 1.
 *
 *   gsl_sf_eta(x) computes the eta function:
 *       eta(x) = (1 - 2^(1-x)) zeta(x).
 *
 *   If optional argument ERR is true, the result, says Y, has an
 *   additional dimension of length 2 prepended to the dimension list of X
 *   which is used to provide an estimate of the error:
 *       Y(1,..) = value of F(X)
 *       Y(2,..) = error estimate for the value of F(X)
 *
 *
 * SEE ALSO: gsl_sf.
 */

extern gsl_poly_solve_quadratic;
extern gsl_poly_solve_cubic;
/* DOCUMENT x = gsl_poly_solve_quadratic(a, b, c);
         or x = gsl_poly_solve_quadratic(v);
         or x = gsl_poly_solve_cubic(a, b, c);
         or x = gsl_poly_solve_cubic(v);

      These functions return the real roots of a quadratic or cubic polynomials
      with real coefficients A, B and C.  When called with a single argument,
      it must be a vector of coefficients: V = [A,B,C].

      If there are no roots, an empty result is returned otherwise a vector of
      1, 2, or 3 roots is returned.  The roots are sorted in ascending order.
      The case of coincident roots is not considered special.  Therefore,
      either 0 or 1 or 2 roots are returned for a quadratic polynomial (a
      single root only occurs if A=0) and either 0 or 1 or 3 roots are
      returned for a cubic polynomial.

      The roots X are such that:

        A*X^2 + B*X + C = 0
        X^3 + A*X^2 + B*X + C = 0

      for a quadratic and a cubic polynomial respectively.

*/

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * End:
 */
