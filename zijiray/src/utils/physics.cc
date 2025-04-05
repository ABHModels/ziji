// clang-format off
#include "utils/physics.h"
#include "utils/novikov_thorne.h"
// clang-format on

namespace ziji {

namespace utils {

namespace physics {

double CutoffPowerlaw(double energy, double gamma, double energy_cut) {
  return pow(energy, -gamma) * exp(-energy / energy_cut);
}

double Powerlaw(double energy, double gamma) { return pow(energy, -gamma); }

double BlackBodyIntensity(double energy_keV, double kbTeff_keV,
                          double f_color_factor) {
  double hplank = phys_constants::h_PLANK_CGS;
  double clight = phys_constants::C_LIGHT_CGS;
  double keV = phys_constants::KEV_erg;

  double A1;
  double Y_angle_func = 1.; // isotropic emission

  A1 = 2.0 / pow(f_color_factor, 4) * pow(keV, 3) / pow(hplank, 3) /
       pow(clight, 2);

  // return I_E
  return A1 * pow(energy_keV, 3.) * Y_angle_func /
         (exp(energy_keV / (f_color_factor * kbTeff_keV)) - 1.0);
}

namespace relativity {

/* get RMS (ISCO) for the Kerr Case */
double KerrRms(double a) {
  //	 accounts for negative spin
  double sign = 1.0;
  if (a < 0) {
    sign = -1.0;
  }

  double Z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) *
                        (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0));
  double Z2 = sqrt((3.0 * a * a) + (Z1 * Z1));

  return 3.0 + Z2 - sign * sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)));
}

double KerrLmom(double r, double a) {
  double r2 = r * r;
  double r0_5 = sqrt(r);
  double r1_5 = r0_5 * r;

  double A = r * r - 2 * a * r0_5 + a * a;
  double B = pow(r, 3. / 4.) * sqrt(r1_5 - 3. * r0_5 + 2 * a);

  return A / B;
}

double KerrEcirc(double r, double a) {
  double u = 1 / r;
  double u2 = u * u;
  double u3 = u2 * u;

  double A = 1. - 2. * u + a * sqrt(u3);
  double B = sqrt(1.0 - 3. * u + 2. * a * sqrt(u3));

  return A / B;
}

double RedshiftLp_Disk(double r, double a, double h) {
  double ut_d =
      ((r * sqrt(r) + a) / (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
  double ut_h = sqrt((h * h + a * a) / (h * h - 2 * h + a * a));

  double gi = ut_d / ut_h;

  return gi;
}

double RedshiftDisk_Screen(double re, double impact_par,
                           const ziji::astro_objects::BlackHole &blackHoleObj) {
  double Omega_var, denom;
  double m[4][4], dmdr[4][4], dmdr2[4][4];

  // Create a Metric object using the BlackHole object
  ziji::GR::Metric metric(blackHoleObj);

  // Compute the metric and its derivatives at r = re, theta = π/2
  metric.computeMetricDerivatives(re, M_PI / 2.0, m, dmdr, dmdr2);

  // Compute angular velocity (Omega_var)
  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];

  // Compute denominator for the redshift calculation
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));

  // Final redshift factor
  double g = denom / (1.0 - impact_par * Omega_var);

  return g;
}

double Loren_dA_dr(double r, double rin,
                   const ziji::astro_objects::BlackHole &blackHoleObj) {
  double Omega, LorentzF, dA_dr;
  double g[4][4], dg_dr[4][4], dg_dr2[4][4];

  // Create a Metric object using the BlackHole object
  ziji::GR::Metric metric(blackHoleObj);

  // Compute the metric and its derivatives at r and theta = π/2
  metric.computeMetricDerivatives(r, M_PI / 2.0, g, dg_dr, dg_dr2);

  if (r < rin) {
    return -1;
  }

  // Compute Omega
  Omega = (-dg_dr[0][3] +
           sqrt(dg_dr[0][3] * dg_dr[0][3] - dg_dr[0][0] * dg_dr[3][3])) /
          dg_dr[3][3];

  // Compute Lorentz factor
  LorentzF = pow(pow(Omega + g[0][3] / g[3][3], 2) * g[3][3] * g[3][3] /
                         (g[0][0] * g[3][3] - g[0][3] * g[0][3]) +
                     1.0,
                 -0.5);

  // Compute the differential area element
  dA_dr = 2 * M_PI * sqrt(g[1][1] * g[3][3]);

  // Return Lorentz factor multiplied by the differential area element
  return LorentzF * dA_dr;
}

double RedshiftDisk_Disk(double re, double ro, double impact_par,
                         const ziji::astro_objects::BlackHole &blackHoleObj) {
  double Omega_var, denom;
  double m[4][4], dmdr[4][4], dmdr2[4][4];
  double kue, kuo;

  // Create a Metric object using the BlackHole object
  ziji::GR::Metric metric(blackHoleObj);

  // Compute the metric and its derivatives at re and theta = π/2
  metric.computeMetricDerivatives(re, M_PI / 2.0, m, dmdr, dmdr2);

  // Compute Omega and denom for the first metric (at re)
  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));

  // Compute kue
  kue = (1 - impact_par * Omega_var) / denom;

  // Compute the metric and its derivatives at ro and theta = π/2
  metric.computeMetricDerivatives(ro, M_PI / 2.0, m, dmdr, dmdr2);

  // Compute Omega and denom for the second metric (at ro)
  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));

  // Compute kuo
  kuo = (1 - impact_par * Omega_var) / denom;

  // Return the ratio of kuo to kue
  return kuo / kue;
}

double
EmissionAngleDisk_Screen(double re, double k_th, double impact_par,
                         const ziji::astro_objects::BlackHole &blackHoleObj) {
  double Omega_var, denom;
  double m[4][4], dmdr[4][4], dmdr2[4][4];

  // Create a Metric object using the BlackHole object
  ziji::GR::Metric metric(blackHoleObj);

  // Compute the metric and its derivatives at r = re, theta = π/2
  metric.computeMetricDerivatives(re, M_PI / 2.0, m, dmdr, dmdr2);

  // Compute angular velocity (Omega_var)
  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];

  // Compute denominator for the redshift calculation
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));

  // redshift factor
  double g = denom / (1.0 - impact_par * Omega_var);

  // cos(theta_emis)
  double cosem = g * sqrt(m[2][2]) * k_th / 1.;

  return fabs(cosem);
  ;
}

} // namespace relativity

} // namespace physics

} // namespace utils

} // namespace ziji
