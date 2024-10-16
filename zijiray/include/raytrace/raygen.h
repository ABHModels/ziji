#ifndef RAYGEN_H
#define RAYGEN_H

// clang-format off
#include "objects/astro_objects.h"
#include "GR/geodesic.h"


namespace ziji {
namespace raytrace {

class RayGen {
public:
  double impact_par;                         // Public attribute for impact parameter
  double vars[NUM_VARS];       // Public attribute for variables (r, theta, etc.)

  // Constructor that takes references to a BlackHole and an AccretionDiskGeometry
  RayGen(const astro_objects::BlackHole &bh, 
         const astro_objects::AccretionDiskGeometry &diskGeom);

  // Public functions to generate rays from different sources
  void FromScreen(double distance, double inclination, double x_scr, double y_scr); // Generate ray from the screen
  void FromDisk(double r_scr, double th_scr, double ph_scr); // Generate ray from the disk
  void FromLampPost(double height, double delta_angle); // Generate ray from the lamp post

private:
  const astro_objects::BlackHole &blackHole_;       // Reference to the BlackHole object
  const astro_objects::AccretionDiskGeometry &diskGeom_; // Reference to the AccretionDiskGeometry
};

} // namespace raytrace
} // namespace ziji

#endif // RAYGEN_H

// clang-format on