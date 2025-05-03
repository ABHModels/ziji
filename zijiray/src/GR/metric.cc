#include "GR/metric.h"

namespace ziji {
namespace GR {

// Constructor implementation
Metric::Metric(const astro_objects::BlackHole &bh) : blackHole(bh) {}

// Compute metric function
void Metric::computeMetric(double r, double th, double (*g)[4]) const {
  double gtt, grr, gthth, gpp, gtp;
  double chi = blackHole.spin; // Use spin from BlackHole object

  double t3 = pow(chi, 2);
  double t2 = pow(r, 2);
  double t5 = cos(th);
  double t6 = pow(t5, 2);
  double t7 = t3 * t6;
  double t8 = t2 + t7;
  double t9 = 1 / t8;
  double t19 = sin(th);
  double t20 = pow(t19, 2);

  gtt = -1 + 2 * r * t9;
  grr = t8 / ((-2 + r) * r + t3);
  gthth = t8;
  gpp = t20 * (t2 + t3 + 2 * r * t20 * t3 * t9);
  gtp = -2 * chi * r * t20 * t9;

  g[0][0] = gtt;
  g[0][3] = gtp;
  g[1][1] = grr;
  g[2][2] = gthth;
  g[3][0] = g[0][3];
  g[3][3] = gpp;
}

// Compute derivatives
void Metric::computeMetricDerivatives(double r, double th, double (*g)[4],
                                      double (*dg)[4], double (*dg2)[4]) const {

  double gtt, grr, gthth, gpp, gtp, dgttdr, dgtpdr, dgppdr, dgttdr2, dgtpdr2,
      dgppdr2;
  double chi = blackHole.spin;

  double t3 = pow(chi, 2);
  double t2 = pow(r, 2);
  double t5 = cos(th);
  double t6 = pow(t5, 2);
  double t7 = t3 * t6;
  double t8 = t2 + t7;
  double t9 = 1 / t8;
  double t19 = sin(th);
  double t20 = pow(t19, 2);
  double t25 = pow(t8, -2);
  double t38 = pow(r, 3);
  double t39 = pow(t8, -3);

  gtt = -1 + 2 * r * t9;
  grr = t8 / ((-2 + r) * r + t3);
  gthth = t8;
  gpp = t20 * (t2 + t3 + 2 * r * t20 * t3 * t9);
  gtp = -2 * chi * r * t20 * t9;

  dgttdr = -4 * t2 * t25 + 2 * t9;
  dgppdr = t20 * (2 * r - 4 * t2 * t20 * t25 * t3 + 2 * t20 * t3 * t9);
  dgtpdr = 4 * chi * t2 * t20 * t25 - 2 * chi * t20 * t9;
  dgttdr2 = -12 * r * t25 + 16 * t38 * t39;
  dgppdr2 = t20 * (2 - 12 * r * t20 * t25 * t3 + 16 * t20 * t3 * t38 * t39);
  dgtpdr2 = 12 * chi * r * t20 * t25 - 16 * chi * t20 * t38 * t39;

  g[0][0] = gtt;
  g[0][3] = gtp;
  g[1][1] = grr;
  g[2][2] = gthth;
  g[3][0] = g[0][3];
  g[3][3] = gpp;

  dg[0][0] = dgttdr;
  dg[0][3] = dgtpdr;
  dg[3][0] = dg[0][3];
  dg[3][3] = dgppdr;

  dg2[0][0] = dgttdr2;
  dg2[0][3] = dgtpdr2;
  dg2[3][0] = dg2[0][3];
  dg2[3][3] = dgppdr2;
}

} // namespace GR
} // namespace ziji