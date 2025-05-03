#ifndef PHYSICS_H
#define PHYSICS_H
#include <cmath>

// clang-format off
#include "utils/phys_constants.h"
#include "GR/metric.h"
// clang-format on

namespace ziji {

namespace utils {

namespace physics {

double CutoffPowerlaw(double energy, double gamma, double energy_cut);
double Powerlaw(double energy, double gamma);
// Cosimo 2012 Eq:3
double BlackBodyIntensity(double energy_keV, double kbTeff_keV,
                          double f_color_factor);

namespace relativity {

double KerrRms(double a);
double KerrLmom(double r, double a);
double KerrEcirc(double r, double a);
double RedshiftLp_Disk(double r, double a, double h);
double RedshiftDisk_Screen(double re, double impact_par,
                           const ziji::astro_objects::BlackHole &blackHoleObj);
double NTFlux_dimless(double re, double rstart,
                      const ziji::GR::Metric &metricObj);
double Loren_dA_dr(double r, double rin,
                   const ziji::astro_objects::BlackHole &blackHoleObj);
double RedshiftDisk_Disk(double re, double ro, double impact_par,
                         const ziji::astro_objects::BlackHole &blackHoleObj);

double
EmissionAngleDisk_Screen(double re, double k_th, double impact_par,
                         const ziji::astro_objects::BlackHole &blackHoleObj);

} // namespace relativity

} // namespace physics

} // namespace utils

} // namespace ziji

#endif // PHYSICS_H