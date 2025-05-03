
// clang-format off
#include "GR/metric.h"
#include "raytrace/raypath.h"

// #include <iostream>
// clang-format on

namespace ziji {
namespace raytrace {

// Initialize static members for default tolerances
double RayPath::defaultAbsTol_[NUM_VARS] = {0, 1.0e-7, 0, 0, 1.0e-7};
double RayPath::defaultRelTol_[NUM_VARS] = {1.0e-7, 0, 1.0e-7, 1.0e-7, 0};

// Constructor: Initializes the RayPath object with default or custom tolerances
RayPath::RayPath(double (&vars)[NUM_VARS], double kphi, double lambdaStep0,
                 const astro_objects::BlackHole &bh,
                 const astro_objects::AccretionDiskGeometry &diskGeom,
                 const double (&absTol)[NUM_VARS],
                 const double (&relTol)[NUM_VARS])
    : vars_(vars), blackHole_(bh), diskGeom_(diskGeom),
      lambdaStep_(lambdaStep0),
      geodesic_(bh, impact_param, absTol,
                relTol) { // Pass absTol, relTol to Geodesic
  // Compute conserved quantities like impact_param
  ComputeConstants_(kphi);
}

// Compute the constants of motion, like impact_param
void RayPath::ComputeConstants_(double kphi) {
  double gm[4][4];                       // Metric components
  double g_tt, g_tp, g_rr, g_thth, g_pp; // Metric components
  double E, kt; // Energy, k_t, k_phi (assuming k_phi = 1.0 for now)

  // Create a temporary Metric object when needed
  ziji::GR::Metric metric(blackHole_);

  // Compute the metric at the initial radius and theta
  metric.computeMetric(vars_[0], vars_[1], gm); // Use Metric class

  // Extract the components of the metric
  g_tt = gm[0][0];
  g_tp = gm[0][3];
  g_rr = gm[1][1];
  g_thth = gm[2][2];
  g_pp = gm[3][3];

  // Compute the energy (E) from the initial variables
  E = sqrt(g_tp * g_tp * kphi * kphi -
           g_tt * (g_rr * vars_[2] * vars_[2] + g_thth * vars_[3] * vars_[3] +
                   g_pp * kphi * kphi));

  // Compute k_t (kt)
  kt = -(g_tp * kphi + E) / g_tt;

  // Compute the impact parameter
  impact_param = -(g_pp * kphi + g_tp * kt) / (g_tt * kt + g_tp * kphi);

  // Normalize initial conditions by energy
  vars_[2] /= E;
  vars_[3] /= E;
}

// Trace the ray and solve geodesic equations
void RayPath::TraceToDisk(double r_stop) {

  double lambda = 0.0; // Initialize lambda (affine parameter)
  int check_disk_under_over;

  if (vars_[1] < M_PI / 2.)
    check_disk_under_over = -1; // under
  else
    check_disk_under_over = 1; // over

  // Solve the geodesic equations using the Geodesic object
  while (fabs(lambdaStep_) > 1.0e-25 &&
         StopAtDisk_(r_stop, check_disk_under_over, lambda)) {
    // std::cout << "r = " << vars_[0] << "\n";
    geodesic_.Solve(vars_, lambda, lambdaStep_);
  }
}

// Private function to stop the geodesic when it hits the disk
bool RayPath::StopAtDisk_(double d, int &check_disk_under_over,
                          double &lambda) {
  bool hor =
      (vars_[0] >
       (2. * defaultRelTol_[0] + 1.) *
           blackHole_.Horizon()); // Horizon check using BlackHole::Horizon()
  bool maxdist = (vars_[0] < d);  // Max distance check
  bool eq_plane = true;           // Theta plane check
  double hcross;
  // Step 1: Check if crossing above or below the disk
  if (check_disk_under_over == -1) {
    if (vars_[1] > 90 * M_PI / 180.) {
      check_disk_under_over = 1;

      // If inside the accretion disk region
      if (vars_[0] < 1.3 * diskGeom_.outer_radius) {
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
      }
      eq_plane = false;
    }
  } else {
    if (vars_[1] < 90 * M_PI / 180.) {
      check_disk_under_over = -1;

      // If inside the accretion disk region
      if (vars_[0] < 1.3 * diskGeom_.outer_radius) {
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
        hcross = -(vars_[1] - M_PI / 2.) / vars_[3];
        geodesic_.Solve(vars_, lambda, hcross);
      }
      eq_plane = false;
    }
  }

  // Step 2: Return result based on horizon, distance, and plane conditions
  return hor * maxdist * eq_plane;
}

} // namespace raytrace
} // namespace ziji