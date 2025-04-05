#ifndef F_ROOT_H
#define F_ROOT_H

#include <cmath>
#include <iostream>

namespace ziji {
namespace utils {
namespace numerical {

double scMethod(double (*f)(double, double *), double *params, double x0,
                double x1, double epsilon, int maxIterations) {
  double x2, f0, f1, f2;
  int iteration = 0;

  f0 = f(x0, params);
  f1 = f(x1, params);

  do {

    /*
    // Check if f(x1) and f(x0) are close to zero
    if (fabs(f1 - f0) < epsilon/1.0e20) {
        cout << "Sc method failed. Division by zero." << endl;
        return NAN;
    }
    */

    x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0);
    f2 = f(x2, params);

    f0 = f1;
    f1 = f2;

    x0 = x1;
    x1 = x2;

    iteration++;
    // cout<<iteration<<": "<<x2<<endl;
  } while (fabs(f2) > epsilon && iteration < maxIterations);

  if (iteration == maxIterations) {
    std::cout << "Sc method failed to converge."
              << "\n";
    return NAN;
  }

  return x2;
}

void scMethod2(double (*f)(double, double *), double *params, double *out,
               double x0, double x1, double epsilon, int maxIterations) {
  double x2, f0, f1, f2;
  int iteration = 0;

  f0 = f(x0, params);
  f1 = f(x1, params);

  // cout<<x0<<"  "<<f0<<endl;
  // cout<<x1<<"  "<<f1<<endl;
  do {

    /*
    // Check if f(x1) and f(x0) are close to zero
    if (fabs(f1 - f0) < epsilon/1.0e20) {
        cout << "Sc method failed. Division by zero." << endl;
    }
    */

    x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0);
    f2 = f(x2, params);

    f0 = f1;
    f1 = f2;

    x0 = x1;
    x1 = x2;

    iteration++;
    // cout<<iteration<<": "<<x2<<"  "<<f2<<endl;
  } while (fabs(f0) > epsilon && iteration < maxIterations);

  if (iteration == maxIterations) {
    std::cout << "Sc method failed to converge."
              << "\n";
  }

  out[0] = x1;                    // x
  out[1] = f1;                    // f
  out[2] = (f1 - f0) / (x1 - x0); //   df/dx = f'(x)
}

double scMethod_singul(double (*f)(double, double *), double *params, double x0,
                       double x1, double epsilon, int maxIterations) {
  double x2, f0, f1, f2;
  int iteration = 0;
  double lambda;

  // Initial function evaluations
  f0 = f(x0, params);
  if (std::isnan(f0)) {

    // std::cout << "Initial point x0 is invalid." << std::endl;
    // return NAN;

    // increase initial point
    lambda = 0.8;
    bool foundValid = false;
    while (lambda > 1e-6) {
      x0 = x0 / lambda; // neeed positive
      f0 = f(x0, params);
      if (!std::isnan(f0)) {
        foundValid = true;
        break;
      }
      lambda /= 2.0;
    }
    if (!foundValid) {
      std::cout << "Unable to find a valid initial point. Exiting."
                << std::endl;
      return NAN;
    }
  }

  f1 = f(x1, params);
  if (std::isnan(f1)) {
    // std::cout << "Initial point x1 is invalid." << std::endl;
    // return NAN;

    // increase initial point
    lambda = 0.5;
    bool foundValid = false;
    while (lambda > 1e-6) {
      x1 = x1 / lambda; // neeed positive
      f1 = f(x1, params);
      if (!std::isnan(f1)) {
        foundValid = true;
        break;
      }
      lambda /= 2.0;
    }
    if (!foundValid) {
      std::cout << "Unable to find a valid initial point. Exiting."
                << std::endl;
      return NAN;
    }
  }

  do {
    // Check if the denominator is too small to avoid division by zero
    if (fabs(f1 - f0) == 0.0) {
      std::cout << "Sc method failed. Division by zero. Iteration: "
                << iteration << std::endl;
      return NAN;
    }

    // Compute the next approximation
    x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0);

    // Evaluate the function at the new point
    f2 = f(x2, params);

    // Check if the new point is valid
    if (std::isnan(f2)) {
      // Reduce step size and retry
      lambda = 0.5;
      bool foundValid = false;
      while (lambda > 1e-6) {
        x2 = x1 + lambda * (x2 - x1);
        f2 = f(x2, params);
        if (!std::isnan(f2)) {
          foundValid = true;
          break;
        }
        lambda /= 2.0;
      }
      if (!foundValid) {
        std::cout << "Unable to find a valid point. Exiting." << std::endl;
        return NAN;
      }
    }

    // Update variables for the next iteration
    f0 = f1;
    f1 = f2;
    x0 = x1;
    x1 = x2;

    iteration++;
    // Uncomment for debugging:
    // std::cout << "Iteration " << iteration << ": x2 = " << x2 << ", f2 = " <<
    // f2
    //           << std::endl;
  } while (fabs(f2) > epsilon && iteration < maxIterations);

  if (iteration == maxIterations) {
    std::cout << "Sc method failed to converge." << std::endl;
    return NAN;
  }

  return x2;
}

} // namespace numerical
} // namespace utils
} // namespace ziji

#endif // F_ROOT_H