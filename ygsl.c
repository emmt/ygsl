/*
 * ygsl.c --
 *
 * Implements Yorick interface to GSL (GNU Scientific Library).
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
 *
 * TO DO:
 *  - gsl_sf_complex_dilog_e fo complex argument
 *  - use gsl_sf_zeta_int, gsl_sf_zetam1_int and
 *    gsl_sf_eta_int for integer argument
 *
 * MISSING:
 *   Dilogarithm with complex argument.
 *   Zeros of Airy Functions.
 *   Zeros of Bessel Functions.
 *   Coulomb Functions.
 *   Legendre Forms of Incomplete Elliptic Integrals.
 *   Carlson Forms of Elliptic Integrals.
 *   Elliptic Functions (Jacobi).
 *   Exponential Functions (some)
 *   Incomplete Fermi-Dirac Integrals
 *   Gamma Function:
 *     - gsl_sf_lngamma_sgn_e, gsl_sf_lngamma_complex_e,
 *     - gsl_sf_fact, gsl_sf_fact_e
 *     - gsl_sf_doublefact, gsl_sf_doublefact_e
 *     - gsl_sf_lnfact, gsl_sf_lnfact_e
 *     - gsl_sf_lndoublefact, gsl_sf_lndoublefact_e
 *     - gsl_sf_choose, gsl_sf_choose_e
 *     - gsl_sf_lnchoose, gsl_sf_lnchoose_e
 *     - gsl_sf_poch, gsl_sf_poch_e
 *     - gsl_sf_lnpoch, gsl_sf_lnpoch_e
 *     - gsl_sf_lnpoch_sgn_e
 *     - gsl_sf_pochrel, gsl_sf_pochrel_e
 *     - gsl_sf_gamma_inc_Q, gsl_sf_gamma_inc_Q_e
 *     - gsl_sf_gamma_inc_P, gsl_sf_gamma_inc_P_e
 *     - gsl_sf_gamma_inc, gsl_sf_gamma_inc_e
 *     - gsl_sf_beta, gsl_sf_beta_e
 *     - gsl_sf_lnbeta, gsl_sf_lnbeta_e
 *     - gsl_sf_beta_inc, gsl_sf_beta_inc_e
 *   Gegenbauer Functions
 *   Hypergeometric Functions
 *   Laguerre Functions
 *   Legendre Functions
 *     - Associated Legendre Polynomials and Spherical Harmonics
 *     - Conical Functions
 *     - Radial Functions for Hyperbolic Space
 *   Logarithmic Functions: gsl_sf_complex_log_e
 *   Power Functions: gsl_sf_pow_int, gsl_sf_pow_int_e
 *   Digamma Functions
 *   Trigonometric Functions: gsl_sf_hypot, ... for complex arguments,
 *       Conversion Functions, Trigonometric Functions With Error Estimates,
 *       Restriction Functions
 *   Zeta Functions:
 *     - gsl_sf_zeta, gsl_sf_zeta_e
 *     - gsl_sf_zetam1_int, gsl_sf_zetam1_int_e
 *     - gsl_sf_hzeta, gsl_sf_hzeta_e
 *     - gsl_sf_eta_int, gsl_sf_eta_int_e
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <yapi.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_poly.h>

/* Define some macros to get rid of some GNU extensions when not compiling
   with GCC. */
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

PLUG_API void y_error(const char *) __attribute__ ((noreturn));

static void fnerr(const char *name, char *msg) __attribute__ ((noreturn));
static void setup(void);
static void error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno);

static int first_time = 1;

static void setup(void)
{
  gsl_set_error_handler(error_handler);
  first_time = 0;
}

static void error_handler(const char *reason,
                          const char *file,
                          int line,
                          int status)
{
  if (reason == NULL || reason[0] == '\0') {
    reason = gsl_strerror(status);
  }
  y_error(reason);
}

/* Same as ypush_d(DIMS) but with an additional dimension of length 2
   prepended to the dimension list DIMS.  Beware that the contents of DIMS
   is modified. */
static double *push_d2(long dims[])
{
  long i;
  if (dims[0] >= (Y_DIMSIZE - 1)) y_error("too many dimensions");
  for (i=dims[0] ; i>=1 ; --i) {
    dims[i + 1] = dims[i];
  }
  dims[1] = 2;
  ++dims[0];
  return ypush_d(dims);
}

static double *get_array_d(int iarg, long *ntot, long dims[])
{
  switch (yarg_typeid(iarg)) {
  case Y_CHAR:
  case Y_SHORT:
  case Y_INT:
  case Y_LONG:
  case Y_FLOAT:
  case Y_DOUBLE:
    return ygeta_d(iarg, ntot, dims);
    break;
  }
  y_error("expecting a non-complex numerical array");
  return NULL; /* avoids compiler warnings */
}

/* 3 - `GSL_PREC_DOUBLE'
     Double-precision, a relative accuracy of approximately 2 * 10^-16.

   2 - `GSL_PREC_SINGLE'
     Single-precision, a relative accuracy of approximately 10^-7.

   1 - `GSL_PREC_APPROX'
     Approximate values, a relative accuracy of approximately 5 * 10^-4.
*/

static void fnerr(const char *name, char *msg)
{
  const size_t bufsiz = 128;
  char buf[bufsiz];

  if (name && name[0]) {
    snprintf(buf, bufsiz, "%s: %s", name, msg);
    buf[bufsiz - 1] = '\0';
    msg = buf;
  }
  y_error(msg);
}

static long get_flags(int iarg, long default_value)
{
  switch (yarg_typeid(iarg)) {
  case Y_VOID: return default_value;
  case Y_CHAR:
  case Y_SHORT:
  case Y_INT:
  case Y_LONG: return ygets_l(iarg);
  }
  y_error("expecting nil or integer scalar");
}

/*---------------------------------------------------------------------------*/

static void sf_driver_1(int argc, const char *name,
                        double (*fn)(double x, gsl_mode_t mode),
                        int (*fn_e)(double x, gsl_mode_t mode,
                                    gsl_sf_result *result))
{
  gsl_sf_result r;
  gsl_mode_t mode;
  long i, ntot, flags, dims[Y_DIMSIZE];
  double *x, *y;

  if (first_time) setup();
  if (argc == 2) {
    flags = get_flags(0, 0L);
    yarg_drop(1);
  } else {
    if (argc != 1) fnerr(name, "takes one or two arguments");
    flags = 0L;
  }
  switch (flags & 0x6) {
  case 2:  mode = GSL_PREC_APPROX; break;
  case 4:  mode = GSL_PREC_SINGLE; break;
  default: mode = GSL_PREC_DOUBLE; break;
  }
  x = get_array_d(0, &ntot, dims);
  if (flags & 0x1) {
    y = push_d2(dims);
    for (i = 0; i < ntot; ++i) {
      fn_e(x[i], mode, &r);
      y[2*i] = r.val;
      y[2*i + 1] = r.err;
    }
  } else {
    y = (yarg_scratch(0) ? x : ypush_d(dims));
    for (i = 0; i < ntot; ++i) {
      y[i] = fn(x[i], mode);
    }
  }
}

#undef FN
#define FN(fn) \
extern void Y_##fn(int argc); \
void Y_##fn(int argc) { sf_driver_1(argc, #fn, fn, fn##_e); }
FN(gsl_sf_airy_Ai)
FN(gsl_sf_airy_Bi)
FN(gsl_sf_airy_Ai_scaled)
FN(gsl_sf_airy_Bi_scaled)
FN(gsl_sf_airy_Ai_deriv)
FN(gsl_sf_airy_Bi_deriv)
FN(gsl_sf_airy_Ai_deriv_scaled)
FN(gsl_sf_airy_Bi_deriv_scaled)
FN(gsl_sf_ellint_Kcomp)
FN(gsl_sf_ellint_Ecomp)
#undef FN

/*---------------------------------------------------------------------------*/

static void sf_driver_2(int argc, const char *name,
                        double (*fn)(double x),
                        int (*fn_e)(double x, gsl_sf_result *result))
{
  gsl_sf_result r;
  long i, ntot, flags, dims[Y_DIMSIZE];
  double *x, *y;

  if (first_time) setup();
  if (argc == 2) {
    flags = get_flags(0, 0L);
    yarg_drop(1);
  } else {
    if (argc != 1) fnerr(name, "takes one or two arguments");
    flags = 0;
  }
  x = get_array_d(0, &ntot, dims);
  if (flags) {
    y = push_d2(dims);
    for (i = 0; i < ntot; ++i) {
      fn_e(x[i], &r);
      y[2*i] = r.val;
      y[2*i+1] = r.err;
    }
  } else {
    y = (yarg_scratch(0) ? x : ypush_d(dims));
    for (i = 0; i < ntot; ++i) {
      y[i] = fn(x[i]);
    }
  }
}

#undef FN
#define FN(fn) \
extern void Y_##fn(int argc); \
void Y_##fn(int argc) { sf_driver_2(argc, #fn, fn, fn##_e); }
FN(gsl_sf_bessel_J0)
FN(gsl_sf_bessel_J1)
FN(gsl_sf_bessel_Y0)
FN(gsl_sf_bessel_Y1)
FN(gsl_sf_bessel_I0)
FN(gsl_sf_bessel_I1)
FN(gsl_sf_bessel_I0_scaled)
FN(gsl_sf_bessel_I1_scaled)
FN(gsl_sf_bessel_K0)
FN(gsl_sf_bessel_K1)
FN(gsl_sf_bessel_K0_scaled)
FN(gsl_sf_bessel_K1_scaled)
FN(gsl_sf_bessel_j0)
FN(gsl_sf_bessel_j1)
FN(gsl_sf_bessel_j2)
FN(gsl_sf_bessel_y0)
FN(gsl_sf_bessel_y1)
FN(gsl_sf_bessel_y2)
FN(gsl_sf_bessel_i0_scaled)
FN(gsl_sf_bessel_i1_scaled)
FN(gsl_sf_bessel_i2_scaled)
FN(gsl_sf_bessel_k0_scaled)
FN(gsl_sf_bessel_k1_scaled)
FN(gsl_sf_bessel_k2_scaled)
FN(gsl_sf_clausen)
FN(gsl_sf_dawson)
FN(gsl_sf_debye_1)
FN(gsl_sf_debye_2)
FN(gsl_sf_debye_3)
FN(gsl_sf_debye_4)
FN(gsl_sf_debye_5)
FN(gsl_sf_debye_6)
FN(gsl_sf_dilog)
FN(gsl_sf_erf)
FN(gsl_sf_erfc)
FN(gsl_sf_log_erfc)
FN(gsl_sf_erf_Z)
FN(gsl_sf_erf_Q)
FN(gsl_sf_hazard)
FN(gsl_sf_exp)
FN(gsl_sf_expm1)
FN(gsl_sf_exprel)
FN(gsl_sf_exprel_2)
FN(gsl_sf_expint_E1)
FN(gsl_sf_expint_E2)
FN(gsl_sf_expint_Ei)
FN(gsl_sf_expint_3)
FN(gsl_sf_Shi)
FN(gsl_sf_Chi)
FN(gsl_sf_Si)
FN(gsl_sf_Ci)
FN(gsl_sf_atanint)
FN(gsl_sf_fermi_dirac_m1)
FN(gsl_sf_fermi_dirac_0)
FN(gsl_sf_fermi_dirac_1)
FN(gsl_sf_fermi_dirac_2)
FN(gsl_sf_fermi_dirac_mhalf)
FN(gsl_sf_fermi_dirac_half)
FN(gsl_sf_fermi_dirac_3half)
FN(gsl_sf_gamma)
FN(gsl_sf_lngamma)
FN(gsl_sf_gammastar)
FN(gsl_sf_gammainv)
FN(gsl_sf_lambert_W0)
FN(gsl_sf_lambert_Wm1)
FN(gsl_sf_legendre_P1)
FN(gsl_sf_legendre_P2)
FN(gsl_sf_legendre_P3)
FN(gsl_sf_legendre_Q0)
FN(gsl_sf_legendre_Q1)
FN(gsl_sf_log)
FN(gsl_sf_log_abs)
FN(gsl_sf_log_1plusx)
FN(gsl_sf_log_1plusx_mx)
FN(gsl_sf_synchrotron_1)
FN(gsl_sf_synchrotron_2)
FN(gsl_sf_transport_2)
FN(gsl_sf_transport_3)
FN(gsl_sf_transport_4)
FN(gsl_sf_transport_5)
FN(gsl_sf_sin)
FN(gsl_sf_cos)
FN(gsl_sf_sinc)
FN(gsl_sf_lnsinh)
FN(gsl_sf_lncosh)
FN(gsl_sf_zeta)
FN(gsl_sf_zetam1)
FN(gsl_sf_eta)
FN(gsl_sf_psi)
FN(gsl_sf_psi_1piy)
FN(gsl_sf_psi_1)
#undef FN

static void sf_driver_3(int argc, const char *name,
                        double (*fn)(int l, double x),
                        int (*fn_e)(int l, double x, gsl_sf_result *result))
{
  gsl_sf_result r;
  long i, ntot, flags, dims[Y_DIMSIZE];
  double *x, *y;
  int l;

  if (first_time) setup();
  if (argc == 3) {
    flags = get_flags(0, 0L);
    yarg_drop(1);
  } else {
    if (argc != 2) fnerr(name, "takes two or three arguments");
    flags = 0;
  }
  l = ygets_l(1);
  x = get_array_d(0, &ntot, dims);
  if (flags) {
    y = push_d2(dims);
    for (i = 0; i < ntot; ++i) {
      fn_e(l, x[i], &r);
      y[2*i] = r.val;
      y[2*i+1] = r.err;
    }
  } else {
    y = (yarg_scratch(0) ? x : ypush_d(dims));
    for (i = 0; i < ntot; ++i) {
      y[i] = fn(l, x[i]);
    }
  }
}

#undef FN
#define FN(fn) \
extern void Y_##fn(int argc); \
void Y_##fn(int argc) { sf_driver_3(argc, #fn, fn, fn##_e); }
FN(gsl_sf_bessel_Jn)
FN(gsl_sf_bessel_Yn)
FN(gsl_sf_bessel_In)
FN(gsl_sf_bessel_In_scaled)
FN(gsl_sf_bessel_Kn)
FN(gsl_sf_bessel_Kn_scaled)
FN(gsl_sf_bessel_jl)
FN(gsl_sf_bessel_yl)
FN(gsl_sf_bessel_il_scaled)
FN(gsl_sf_bessel_kl_scaled)
FN(gsl_sf_exprel_n)
FN(gsl_sf_fermi_dirac_int)
FN(gsl_sf_taylorcoeff)
FN(gsl_sf_legendre_Pl)
FN(gsl_sf_legendre_Ql)
FN(gsl_sf_psi_n)
#undef FN

static void sf_driver_4(int argc, const char *name,
                        double (*fn)(double nu, double x),
                        int (*fn_e)(double nu, double x,
                                    gsl_sf_result *result))
{
  double nu;
  gsl_sf_result r;
  long i, ntot, flags, dims[Y_DIMSIZE];
  double *x, *y;

  if (first_time) setup();
  if (argc == 3) {
    flags = get_flags(0, 0L);
    yarg_drop(1);
  } else {
    if (argc != 2) fnerr(name, "takes one or two arguments");
    flags = 0;
  }
  nu = ygets_d(1);
  x = get_array_d(0, &ntot, dims);
  if (flags) {
    y = push_d2(dims);
    for (i = 0; i < ntot; ++i) {
      fn_e(nu, x[i], &r);
      y[2*i] = r.val;
      y[2*i+1] = r.err;
    }
  } else {
    y = (yarg_scratch(0) ? x : ypush_d(dims));
    for (i = 0; i < ntot; ++i) {
      y[i] = fn(nu, x[i]);
    }
  }
}

#undef FN
#define FN(fn) \
extern void Y_##fn(int argc); \
void Y_##fn(int argc) { sf_driver_4(argc, #fn, fn, fn##_e); }
FN(gsl_sf_bessel_Jnu)
FN(gsl_sf_bessel_Ynu)
FN(gsl_sf_bessel_Inu)
FN(gsl_sf_bessel_Inu_scaled)
FN(gsl_sf_bessel_Knu)
FN(gsl_sf_bessel_lnKnu)
FN(gsl_sf_bessel_Knu_scaled)
#undef FN

/*---------------------------------------------------------------------------*/
/* POLYNOMIAL ROOTS */

static void push_vector_d(long n, const double inp[])
{
  long dims[2];
  double* out;
  long j;

  if (n > 0) {
    dims[0] = 1;
    dims[1] = n;
    out = ypush_d(dims);
    for (j = 0; j < n; ++j) {
      out[j] = inp[j];
    }
  } else {
    ypush_nil();
  }
}

void Y_gsl_poly_solve_quadratic(int argc)
{
  double a, b, c;
  long dims[Y_DIMSIZE];
  double x[2];
  const double* coef;
  long n;

  if (argc == 1) {
    coef = ygeta_d(0, &n, dims);
    if (dims[0] == 1 && dims[1] == 3) {
      a = coef[0];
      b = coef[1];
      c = coef[2];
    } else {
      goto bad_args;
    }
  } else if (argc == 3) {
    a = ygets_d(2);
    b = ygets_d(1);
    c = ygets_d(0);
  } else {
    bad_args:
    y_error("expecting a 3-element vector or 3 arguments");
    return;
  }
  n = gsl_poly_solve_quadratic(a, b, c, &x[0], &x[1]);
  push_vector_d(n, x);
}

void Y_gsl_poly_solve_cubic(int argc)
{
  double a, b, c;
  long dims[Y_DIMSIZE];
  double x[3];
  const double* coef;
  long n;

  if (argc == 1) {
    coef = ygeta_d(0, &n, dims);
    if (dims[0] == 1 && dims[1] == 3) {
      a = coef[0];
      b = coef[1];
      c = coef[2];
    } else {
      goto bad_args;
    }
  } else if (argc == 3) {
    a = ygets_d(2);
    b = ygets_d(1);
    c = ygets_d(0);
  } else {
    bad_args:
    y_error("expecting a 3-element vector or 3 arguments");
    return;
  }
  n = gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
  push_vector_d(n, x);
}

/*
 * Local Variables:
 * mode: C
 * c-basic-offset: 2
 * tab-width: 8
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
