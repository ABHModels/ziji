#ifndef RAYPATH_H
#define RAYPATH_H

// clang-format off

#include "GR/geodesic.h"
// clang-format on

namespace ziji {
namespace raytrace {

class RayPath {
public:
  double impact_param;

  // Constructor with default tolerances for absTol and relTol
  RayPath(double (&vars)[NUM_VARS], double kphi, double lambdaStep0,
          const astro_objects::BlackHole &bh,
          const astro_objects::AccretionDiskGeometry &diskGeom,
          const double (&absTol)[NUM_VARS] = defaultAbsTol_,
          const double (&relTol)[NUM_VARS] = defaultRelTol_);

  // Public function to trace the path to the disk
  void TraceToDisk(double r_stop);

private:
  double (&vars_)[NUM_VARS]; // Reference to vars[] array
  ziji::GR::Geodesic
      geodesic_; // Geodesic object for solving geodesic equations
  const astro_objects::BlackHole
      &blackHole_; // Reference to the BlackHole object
  const astro_objects::AccretionDiskGeometry
      &diskGeom_;     // Accretion Disk Geometry
  double lambdaStep_; // Step size for affine parameter

  // Default tolerances (set as static members)
  static double defaultAbsTol_[NUM_VARS];
  static double defaultRelTol_[NUM_VARS];

  // Compute constants like impact_param and conserved quantities
  void ComputeConstants_(double kphi);
  // Private function that stops the geodesic when it hits the disk
  bool StopAtDisk_(double d, int &check_disk_under_over, double &lambda);
};

} // namespace raytrace
} // namespace ziji

#endif // RAYPATH_H