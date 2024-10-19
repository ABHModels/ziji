// clang-format off
#include "GR/metric.h"
#include "raytrace/raypath.h"
#include "raytrace/raygen.h"
// clang-format on

namespace ziji {
namespace raytrace {

// Constructor initializes the RayGen object
RayGen::RayGen(const astro_objects::BlackHole &bh,
               const astro_objects::AccretionDiskGeometry &diskGeom)
    : blackHole_(bh), diskGeom_(diskGeom) {}

// Generate ray from the screen
void RayGen::FromScreen(double distance, double inclination, double x_scr,
                        double y_scr) {
  double d, dcs;
  double d2, r2, alpha, alpha2, beta, beta2, betas, dcsbetas;
  double kphi;

  d = distance;
  d2 = d * d;
  dcs = d * cos(inclination);
  alpha = x_scr;
  beta = y_scr;

  alpha2 = alpha * alpha;

  beta2 = beta * beta;
  betas = beta * sin(inclination);
  dcsbetas = dcs + betas;

  vars[0] = sqrt(d2 + alpha2 + beta2);
  vars[1] = acos(dcsbetas / vars[0]);

  r2 = vars[0] * vars[0];

  vars[2] = d / vars[0];
  vars[3] =
      (-cos(inclination) + d / r2 * dcsbetas) / sqrt(r2 - dcsbetas * dcsbetas);

  vars[4] = atan2(alpha, (d * sin(inclination) - beta * cos(inclination)));
  kphi = -alpha * sin(inclination) /
         (pow(d * sin(inclination) - beta * cos(inclination), 2.0) + alpha2);

  // Example of a dummy value for lambdaStep0
  double lambdaStep0 = -128.;

  // Creating an instance of RayPath
  RayPath rayPath(vars, kphi, lambdaStep0, blackHole_, diskGeom_);

  // Call the RayPath method, for example, TraceToDisk
  rayPath.TraceToDisk(2 * distance);

  impact_par = rayPath.impact_param;
}

// Generate ray from the disk
void RayGen::FromDisk(double r_scr, double th_scr, double ph_scr) {
  // Set the initial conditions based on the disk position

  double kt, kphi;
  double E;
  double gm[4][4], dgm[4][4], d2gm[4][4];

  double pt, pr, pth, pphi;
  double pt_NR, pr_NR, pth_NR, pphi_NR;
  double e2nu, e2mu, e2sigma, e2lambda, w2, Omega, uet, vm, Gamma, unorm, w;
  double enu_, emu_, elambda_, esigma_;

  // Coord transform ////////////////////////

  pt = 1.;
  pr = -sin(th_scr) * cos(ph_scr);
  pth = cos(th_scr);
  pphi = -sin(th_scr) * sin(ph_scr);

  ziji::GR::Metric metric(blackHole_);
  metric.computeMetricDerivatives(r_scr, M_PI / 2, gm, dgm, d2gm);

  e2nu = -gm[0][0] + gm[0][3] * gm[0][3] / gm[3][3];
  e2mu = gm[1][1];
  e2lambda = gm[2][2];
  e2sigma = gm[3][3];
  w = -gm[0][3] / gm[3][3];
  enu_ = 1.0 / sqrt(e2nu);
  emu_ = 1.0 / sqrt(e2mu);
  elambda_ = 1.0 / sqrt(e2lambda);
  esigma_ = 1.0 / sqrt(e2sigma);

  Omega = (-dgm[0][3] + sqrt(dgm[0][3] * dgm[0][3] - dgm[0][0] * dgm[3][3])) /
          dgm[3][3]; // angular velocity

  vm = (Omega - w) * sqrt(e2sigma / e2nu);
  Gamma = 1 / sqrt(1 - vm * vm);

  pt_NR = Gamma * (pt + vm * pphi);
  pr_NR = pr;
  pth_NR = pth;
  pphi_NR = Gamma * (vm * pt + pphi);

  ////////////////

  vars[0] = r_scr;               // r
  vars[1] = M_PI / 2. - 1.0e-10; // th
  vars[4] = 0.;                  // ph

  vars[2] = emu_ * pr_NR;      // kr
  vars[3] = elambda_ * pth_NR; // kth
  kphi = w * enu_ * pt_NR + esigma_ * pphi_NR;

  double lambdaStep0 = -1.0;

  // Creating an instance of RayPath
  RayPath rayPath(vars, kphi, lambdaStep0, blackHole_, diskGeom_);

  // Call the RayPath method, for example, TraceToDisk
  rayPath.TraceToDisk(2 * diskGeom_.outer_radius);

  impact_par = rayPath.impact_param;
}

// Generate ray from the lamp post
void RayGen::FromLampPost(double height, double delta_angle) {
  // Set the initial conditions based on the height of the lamp post and angular
  double kphi;
  double gm[4][4];
  double g_tt, g_tp, g_rr, g_thth, g_pp;
  double unorm;

  vars[0] = height; // r
  vars[1] = 1.0e-8; // th
  vars[4] = 0.1;    // ph

  ziji::GR::Metric metric(blackHole_);
  metric.computeMetric(height, vars[1], gm);

  g_rr = gm[1][1];
  g_thth = gm[2][2];

  vars[2] = -cos(delta_angle) / sqrt(g_rr);  // kr
  vars[3] = sin(delta_angle) / sqrt(g_thth); // kth

  kphi = 0.;

  // Example of a dummy value for lambdaStep0
  double lambdaStep0 = 1.0;

  // Creating an instance of RayPath
  RayPath rayPath(vars, kphi, lambdaStep0, blackHole_, diskGeom_);

  // Call the RayPath method, for example, TraceToDisk
  rayPath.TraceToDisk(2 * diskGeom_.outer_radius);

  impact_par = rayPath.impact_param;
}

} // namespace raytrace
} // namespace ziji