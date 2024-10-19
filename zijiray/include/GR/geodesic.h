#ifndef GR_GEODESIC_H
#define GR_GEODESIC_H

// clang-format off
#include "objects/astro_objects.h"
#include "utils/RKF45.h"  // Assuming you have an implementation of RKF45





namespace ziji {
// Global constant for number of variables
static const int NUM_VARS = 5; // Change to 4 if you want to ignore phi

namespace GR {

class Geodesic {

public:
  double prev_var[NUM_VARS]; // Public variable to store the previous vars

private:
  const astro_objects::BlackHole &blackHole;  // Reference to the BlackHole object
  const double &impact_param;                 // Impact parameter (reference)
  const double (&absTol)[NUM_VARS];           // Absolute tolerance array (reference, size NUM_VARS)
  const double (&relTol)[NUM_VARS];           // Relative tolerance array (reference, size NUM_VARS)

  // RKF45 object for solving ODEs
  utils::numerical::RKF45<NUM_VARS> rkF;

  // Static geodesic equation solver that will accept the class instance as param
  static void GeoEquations_(double vars[], double diffs[], double lambda, double* param);


public:
  // Constructor
  Geodesic(const astro_objects::BlackHole &bh, const double &impact_param,
           const double (&absTol)[NUM_VARS], const double (&relTol)[NUM_VARS]);

  // // Geodesic equation
  // void GeoEquations(double vars[], double diffs[], double lambda,
  //                          double *param);

  // Solve
  void Solve(double vars[NUM_VARS], double &lambda, double &lambdaStep);
};

} // namespace GR
} // namespace ziji

#endif // GR_GEODESIC_H