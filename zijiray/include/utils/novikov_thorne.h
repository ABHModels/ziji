
// clang-format off
#include <gsl/gsl_integration.h>
#include "GR/metric.h"
// clang-format on

// Global static reference (only visible in this translation unit)
static const ziji::GR::Metric *metricObjRef = nullptr;

static double IntegFlux(double re, double rstart);

namespace ziji {
namespace utils {
namespace physics {
namespace relativity {

double NTFlux_dimless(double re, double rstart,
                      const ziji::GR::Metric &metricObj) {
  double E_var, Lz_var, Omega_var, denom, dOmega_var_dr;
  double m[4][4], dmdr[4][4], dmdr2[4][4];

  metricObjRef = &metricObj;
  metricObj.computeMetricDerivatives(re, M_PI / 2.0, m, dmdr, dmdr2);

  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));
  E_var = -(m[0][0] + m[0][3] * Omega_var) / denom;
  Lz_var = (m[0][3] + m[3][3] * Omega_var) / denom;

  dOmega_var_dr =
      -Omega_var * dmdr2[3][3] / dmdr[3][3] +
      (-dmdr2[0][3] +
       (2 * dmdr[0][3] * dmdr2[0][3] - dmdr2[0][0] * dmdr[3][3] -
        dmdr[0][0] * dmdr2[3][3]) /
           2 / sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
          dmdr[3][3];

  double integ, sqdetg; //, g;

  // g = denom/(1.0 - impact_par*Omega_var);
  //*gredsh = g;

  integ = IntegFlux(re, rstart);
  sqdetg = sqrt(m[1][1] * (m[0][3] * m[0][3] - m[0][0] * m[3][3]));

  return -1.0 / sqdetg * dOmega_var_dr / pow(E_var - Omega_var * Lz_var, 2) *
         integ; //*pow(g,4)
}

} // namespace relativity

} // namespace physics

} // namespace utils

} // namespace ziji

static double IntegFunc(double r, void *params);

static double IntegFlux(double re, double rstart) {
  double result, error;
  double p1, p2;

  p1 = rstart;
  p2 = re;

  double alpha = 1.0;
  // gsl_set_error_handler_off();
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  gsl_function F;
  F.function = &IntegFunc;
  F.params = &alpha;

  gsl_integration_qags(&F, p1, p2, 0., 1.0e-7, 1000, w, &result, &error);

  gsl_integration_workspace_free(w); // 10.03.2022

  return result;
}

static double IntegFunc(double r, void *params) {
  double E_var, Lz_var, Omega_var, denom, dLz_var_dr, dOmega_var_dr;
  double m[4][4], dmdr[4][4], dmdr2[4][4];

  metricObjRef->computeMetricDerivatives(r, M_PI / 2., m, dmdr, dmdr2);

  Omega_var =
      (-dmdr[0][3] + sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
      dmdr[3][3];
  denom = sqrt(
      -(m[0][0] + 2.0 * m[0][3] * Omega_var + m[3][3] * Omega_var * Omega_var));
  E_var = -(m[0][0] + m[0][3] * Omega_var) / denom;
  Lz_var = (m[0][3] + m[3][3] * Omega_var) / denom;

  dOmega_var_dr =
      -Omega_var * dmdr2[3][3] / dmdr[3][3] +
      (-dmdr2[0][3] +
       (2 * dmdr[0][3] * dmdr2[0][3] - dmdr2[0][0] * dmdr[3][3] -
        dmdr[0][0] * dmdr2[3][3]) /
           2 / sqrt(dmdr[0][3] * dmdr[0][3] - dmdr[0][0] * dmdr[3][3])) /
          dmdr[3][3];

  double exp1, exp2;
  exp1 = -dmdr[0][0] - 2 * (dmdr[0][3] * Omega_var + m[0][3] * dOmega_var_dr) -
         (dmdr[3][3] * Omega_var * Omega_var +
          m[3][3] * 2 * Omega_var * dOmega_var_dr);
  exp2 = dmdr[0][3] + dmdr[3][3] * Omega_var + m[3][3] * dOmega_var_dr;

  // True version 25.02.2022
  dLz_var_dr = -Lz_var / (2 * denom * denom) * exp1 + exp2 / denom;
  //////////////////////////

  return (E_var - Omega_var * Lz_var) * dLz_var_dr;
}
