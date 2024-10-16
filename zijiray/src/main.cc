

#include "objects/astro_objects.h"
#include "raytrace/raygen.h"
#include <iostream>

int main() {
  // Create a BlackHole object
  ziji::astro_objects::BlackHole bh(0.8, 10.); // Example spin value

  // Create an AccretionDiskGeometry object
  ziji::astro_objects::AccretionDiskGeometry diskGeom(
      1.0,
      500.0); // Example disk radii

  // Initialize the RayGen object
  ziji::raytrace::RayGen rayGen(bh, diskGeom);

  // Example: Generating a ray from the screen
  double distance = 10000.0;             // Example distance
  double inclination = M_PI * 60. / 180; // Example inclination (in radians)
  double x_scr = 5.;                     // Screen coordinates (example)
  double y_scr = 0.0001;                 // Screen coordinates (example)

  // Generate ray from screen
  rayGen.FromScreen(distance, inclination, x_scr, y_scr);

  std::cout << "rem = " << rayGen.vars[0] << "\n";

  return 0;
}