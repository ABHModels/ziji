#ifndef ASTRO_OBJECTS_H
#define ASTRO_OBJECTS_H

// #include <Eigen/Dense>
#include <cmath>
#include <optional>
#include <vector>

namespace ziji {

namespace astro_objects {

// using ArrayXd = Eigen::ArrayXd; // Alias for Eigen's array type

struct BlackHole {
public:
  double spin; // Dimensionless spin
  double mass; // Mass in solar masses

  // Constructor
  BlackHole(double spin_, double mass_) : spin(spin_), mass(mass_), rhor_(1.) {
    calculHorizon(); // Calculate the horizon during construction
  }

  // Initialize the black hole with specific spin and mass
  void init(double spin_, double mass_) {
    spin = spin_;
    mass = mass_;
    calculHorizon();
  }

  // Public function to return the horizon radius
  double Horizon() const {
    return rhor_; // Return the previously calculated horizon radius
  }

private:
  double rhor_; // Private variable to store the horizon radius

  // Private function to calculate the horizon radius
  void calculHorizon() {
    // Example formula for calculating the horizon radius
    rhor_ = 1. + sqrt(1. - spin * spin); // Simplified Kerr radius
  }
};

// Struct for the Corona with optional scale_luminosity
struct Corona {
  double gamma;  // Photon index of the corona
  double ecut;   // Cutoff energy (keV)
  double height; // Height in gravitational radii (rg)
  std::optional<double>
      scale_luminosity; // in L_edington, dimless, for example 0.1 * Led

  // Constructor with optional scale_luminosity
  Corona(double gamma, double ecut, double height,
         std::optional<double> scale_luminosity = std::nullopt)
      : gamma(gamma), ecut(ecut), height(height),
        scale_luminosity(scale_luminosity) {}
};

// Struct for Accretion Disk Geometry (size)
struct AccretionDiskGeometry {
  double inner_radius; // Inner radius of the disk
  double outer_radius; // Outer radius of the disk

  // Constructor
  AccretionDiskGeometry(double inner_r, double outer_r)
      : inner_radius(inner_r), outer_radius(outer_r) {}
};

// Struct for Accretion Disk Properties (physical parameters)
// struct AccretionDiskProperties {
//   // Optional parameters
//   std::optional<double> density;          // Optional density
//   std::optional<double> ionization;       // Optional ionization
//   std::optional<double> scale_luminosity; // in L_edington, dimensionless
//   std::optional<double> iron_abundance;   // Optional iron abundance
//   std::optional<double> color_factor;     // Optional color factor

//   std::vector<int> profile = {
//       0, 0}; // Array to indicate constant density/ionization

//   // Optional arrays for density and ionization profiles
//   std::optional<ArrayXd> density_profile;    // Optional density profile
//   std::optional<ArrayXd> ionization_profile; // Optional ionization profile

//   // Constructor
//   AccretionDiskProperties(
//       std::vector<int> profile = {0, 0},
//       std::optional<double> density = std::nullopt,
//       std::optional<double> ionization = std::nullopt,
//       std::optional<double> scale_luminosity = std::nullopt,
//       std::optional<double> iron_abundance = std::nullopt,
//       std::optional<double> color_factor = std::nullopt,
//       std::optional<ArrayXd> density_profile = std::nullopt,
//       std::optional<ArrayXd> ionization_profile = std::nullopt)
//       : profile(profile), density(density), ionization(ionization),
//         scale_luminosity(scale_luminosity), iron_abundance(iron_abundance),
//         color_factor(color_factor), density_profile(density_profile),
//         ionization_profile(ionization_profile) {}
// };

} // namespace astro_objects
} // namespace ziji

#endif // ASTRO_OBJECTS_H