#ifndef F_ROOT_H
#define F_ROOT_H

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

} // namespace numerical
} // namespace utils
} // namespace ziji

#endif // F_ROOT_H