#include "GR/geodesic.h"
#include <iostream>

namespace ziji {
namespace GR {

// Constructor that initializes everything, including rkF
Geodesic::Geodesic(const astro_objects::BlackHole &bh,
                   const double &impact_param, const double (&absTol)[NUM_VARS],
                   const double (&relTol)[NUM_VARS])
    : blackHole(bh), impact_param(impact_param), absTol(absTol), relTol(relTol),
      rkF(absTol, relTol, GeoEquations_, reinterpret_cast<double *>(this)) {
  // Initialize prev_var with zeros or initial values
  for (int i = 0; i < NUM_VARS; ++i) {
    prev_var[i] = 0.0;
  }
}

// Solve the geodesic equations using RKF45
void Geodesic::Solve(double vars[NUM_VARS], double &lambda,
                     double &lambdaStep) {
  rkF.OdeSolve(prev_var, vars, &lambda, &lambdaStep);
}

// Geodesic equations
void Geodesic::GeoEquations_(double vars[], double diffs[], double lambda,
                             double *param) {
  double r, th;
  double kt, kphi;
  double kt2, kr2, kth2, kp2, ktp, krth;
  double r2, r3, r4, r5, r6, chi2, chi3, chi4, cs, cs2, cs4, s, s2, css;
  double delta, rho2, rho4, rho6, rho2neg, sum, tr, chi2css, chi2s2, chis2;
  double denom1, denom2, deltarho2neg;
  double ch_rtt, ch_rtp, ch_rrr, ch_rrth, ch_rthth, ch_rpp;
  double ch_thtt, ch_thtp, ch_thrr, ch_thrth, ch_ththth, ch_thpp;
  double g_tt, g_tp, g_pp;

  // Cast param back to Geodesic instance
  Geodesic *instance = reinterpret_cast<Geodesic *>(param);

  double impact_par = instance->impact_param;

  double chi = instance->blackHole.spin;

  r = vars[0];
  th = vars[1];

  r2 = r * r;
  r3 = r2 * r;
  r4 = r3 * r;
  r5 = r4 * r;
  r6 = r5 * r;

  tr = 2.0 * r;

  chi2 = chi * chi;
  chi3 = chi2 * chi;
  chi4 = chi3 * chi;

  cs = cos(th);
  cs2 = cs * cs;
  cs4 = cs2 * cs2;
  s = sin(th);
  s2 = s * s;
  css = cs * s;

  chis2 = chi * s2;
  chi2css = chi2 * css;
  chi2s2 = chi2 * s2;

  sum = r2 + chi2;

  delta = sum - tr;

  rho2 = r2 + chi2 * cs2;

  rho4 = rho2 * rho2;

  rho6 = rho4 * rho2;

  rho2neg = rho2 - 2.0 * r2;

  deltarho2neg = delta * rho2neg;

  denom1 = rho2 * delta;

  ch_rtt = -deltarho2neg / rho6;

  ch_rtp = -ch_rtt * chis2;

  ch_rrr = (chi2 * (cs2 + r * s2) - r2) / denom1;

  ch_rrth = -chi2css / rho2;

  ch_rthth = -r * delta / rho2;

  ch_rpp =
      delta * s2 *
      (chi4 * (1.0 - r) * cs4 - chi2 * (cs2 * (1.0 + tr * r2) - r2 * s2) - r5) /
      rho6;

  ch_thtt = -tr * chi2css / rho6;

  ch_thtp = tr * chi * css * sum / rho6;

  ch_thrr = chi2css / denom1;

  ch_thrth = r / rho2;

  ch_ththth = -chi2css / rho2;

  ch_thpp = -css *
            (chi4 * delta * cs4 + 2.0 * chi2 * r2 * delta * cs2 + tr * chi4 +
             (4.0 * r3 + r4) * chi2 + r6) /
            rho6;

  g_tt = -1.0 + tr / rho2;

  g_tp = -tr * chi * s2 / rho2;

  g_pp = s2 * (sum + tr * chi2s2 / rho2);

  denom2 = (g_tt * g_pp - g_tp * g_tp);

  kt = -(g_pp + impact_par * g_tp) / denom2;

  kphi = (g_tp + impact_par * g_tt) / denom2;

  diffs[0] = vars[2];

  diffs[1] = vars[3];

  kt2 = kt * kt;
  kr2 = vars[2] * vars[2];
  kth2 = vars[3] * vars[3];
  kp2 = kphi * kphi;
  ktp = kt * kphi;
  krth = vars[2] * vars[3];

  diffs[2] = -(ch_rtt * kt2 + ch_rrr * kr2 + ch_rthth * kth2 + ch_rpp * kp2 +
               2.0 * (ch_rtp * ktp + ch_rrth * krth));

  diffs[3] = -(ch_thtt * kt2 + ch_thrr * kr2 + ch_ththth * kth2 +
               ch_thpp * kp2 + 2.0 * (ch_thtp * ktp + ch_thrth * krth));

  // if (NUM_VARS == 5) {
  diffs[4] = kphi; // dphi/dlambda (only if NUM_VARS == 5)
  // }
}

} // namespace GR
} // namespace ziji