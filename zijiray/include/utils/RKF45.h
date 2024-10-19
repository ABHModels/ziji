#ifndef RKF45_H
#define RKF45_H

#include <cmath>

namespace ziji {
namespace utils {
namespace numerical {

template <int Nmax_ODE_RK> class RKF45 {

public:
  RKF45(const double (&absolutTol)[Nmax_ODE_RK],
        const double (&relativeTol)[Nmax_ODE_RK],
        void func_eq(double *, double *, double, double *), double *par_eq);

  // Solve ODE
  void OdeSolve(double y0_[], double y_[], double *x, double *h);

  void rkStepper(double y_[], double dydx0_[], double *x, double *h);

private:
  // Private variables
  const double (&absTol_)[Nmax_ODE_RK]; // Absolute error tolerance (reference)
  const double (&relTol_)[Nmax_ODE_RK]; // Relative error tolerance (reference)

  int N_ode; // Number of ODEs

  // RKF45 Constants
  static const double c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51,
      a52, a53, a54, a61, a62, a63, a64, a65, b1, b3, b4, b5, b6, e1, e3, e4,
      e5, e6;

  // Function pointer declaration
  void (*equation_ode)(double *, double *, double,
                       double *); //  (*f)(y[], dydx[], x, *params)
  double *params_eq;              // Parameters for equation_ode

  double normaliseError(double y0[], double y[], double err[]);

  double maxAB(double A, double B);
};

// clang-format off
// Definition of static constants outside the class
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::c2 = 1.0 / 4.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::c3 = 3.0 / 8.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::c4 = 12.0 / 13.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::c5 = 1.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::c6 = 1.0 / 2.0;

template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a21 = 1.0 / 4.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a31 = 3.0 / 32.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a32 = 9.0 / 32.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a41 = 1932.0 / 2197.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a42 = -7200.0 / 2197.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a43 = 7296.0 / 2197.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a51 = 439.0 / 216.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a52 = -8.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a53 = 3680.0 / 513.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a54 = -845.0 / 4104.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a61 = -8.0 / 27.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a62 = 2.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a63 = -3544.0 / 2565.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a64 = 1859.0 / 4104.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::a65 = -11.0 / 40.0;

template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::b1 = 16.0 / 135.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::b3 = 6656.0 / 12825.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::b4 = 28561.0 / 56430.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::b5 = -9.0 / 50.0;
template <int Nmax_ODE_RK> const double RKF45<Nmax_ODE_RK>::b6 = 2.0 / 55.0;

template <int Nmax_ODE_RK>
const double RKF45<Nmax_ODE_RK>::e1 = 25.0 / 216.0 - RKF45<Nmax_ODE_RK>::b1;
template <int Nmax_ODE_RK>
const double RKF45<Nmax_ODE_RK>::e3 = 1408.0 / 2565.0 - RKF45<Nmax_ODE_RK>::b3;
template <int Nmax_ODE_RK>
const double RKF45<Nmax_ODE_RK>::e4 = 2197.0 / 4104.0 - RKF45<Nmax_ODE_RK>::b4;
template <int Nmax_ODE_RK>
const double RKF45<Nmax_ODE_RK>::e5 = -1.0 / 5.0 - RKF45<Nmax_ODE_RK>::b5;
template <int Nmax_ODE_RK>
const double RKF45<Nmax_ODE_RK>::e6 = -RKF45<Nmax_ODE_RK>::b6;
// clang-format on

// Constructor
template <int Nmax_ODE_RK>
RKF45<Nmax_ODE_RK>::RKF45(const double (&absolutTol)[Nmax_ODE_RK],
                          const double (&relativeTol)[Nmax_ODE_RK],
                          void func_eq(double *, double *, double, double *),
                          double *par_eq)
    : absTol_(absolutTol), relTol_(relativeTol), N_ode(Nmax_ODE_RK),
      equation_ode(func_eq), params_eq(par_eq) {}

// Solve ODE
template <int Nmax_ODE_RK>
void RKF45<Nmax_ODE_RK>::OdeSolve(double y0_[], double y_[], double *x,
                                  double *h) {
  double dydx0_[Nmax_ODE_RK];

  for (int i = 0; i < N_ode; i++) {
    y0_[i] = y_[i];
  }

  equation_ode(y_, dydx0_, *x, params_eq); // initial dydx
  rkStepper(y_, dydx0_, x, h);
}

template <int Nmax_ODE_RK>
void RKF45<Nmax_ODE_RK>::rkStepper(double y_[], double dydx0_[], double *x,
                                   double *h) {
  double k1_[Nmax_ODE_RK], k2_[Nmax_ODE_RK], k3_[Nmax_ODE_RK], k4_[Nmax_ODE_RK],
      k5_[Nmax_ODE_RK], k6_[Nmax_ODE_RK], yTemp_[Nmax_ODE_RK],
      dydxTemp_[Nmax_ODE_RK], err_[Nmax_ODE_RK];

  int i;
  double x0 = *x;
  double dx = *h;
  double ep_ep0;

  // Step calculations
  // 1
  for (i = 0; i < N_ode; i++) {
    k1_[i] = dx * dydx0_[i];
    yTemp_[i] = y_[i] + a21 * k1_[i];
  }

  // 2
  equation_ode(yTemp_, dydxTemp_, x0 + c2 * dx, params_eq);
  for (i = 0; i < N_ode; i++) {
    k2_[i] = dx * dydxTemp_[i];
    yTemp_[i] = y_[i] + a31 * k1_[i] + a32 * k2_[i];
  }

  // 3
  equation_ode(yTemp_, dydxTemp_, x0 + c3 * dx, params_eq);
  for (i = 0; i < N_ode; i++) {
    k3_[i] = dx * dydxTemp_[i];
    yTemp_[i] = y_[i] + a41 * k1_[i] + a42 * k2_[i] + a43 * k3_[i];
  }

  // 4
  equation_ode(yTemp_, dydxTemp_, x0 + c4 * dx, params_eq);
  for (i = 0; i < N_ode; i++) {
    k4_[i] = dx * dydxTemp_[i];
    yTemp_[i] =
        y_[i] + a51 * k1_[i] + a52 * k2_[i] + a53 * k3_[i] + a54 * k4_[i];
  }

  // 5
  equation_ode(yTemp_, dydxTemp_, x0 + c5 * dx, params_eq);
  for (i = 0; i < N_ode; i++) {
    k5_[i] = dx * dydxTemp_[i];
    yTemp_[i] = y_[i] + a61 * k1_[i] + a62 * k2_[i] + a63 * k3_[i] +
                a64 * k4_[i] + a65 * k5_[i];
  }

  // 6
  equation_ode(yTemp_, dydxTemp_, x0 + c6 * dx, params_eq);
  for (i = 0; i < N_ode; i++) {
    k6_[i] = dx * dydxTemp_[i];

    // Calculate the 5th-order solution
    yTemp_[i] = y_[i] + b1 * k1_[i] + b3 * k3_[i] + b4 * k4_[i] + b5 * k5_[i] +
                b6 * k6_[i];
  }

  // Error estimation
  for (i = 0; i < N_ode; i++) {
    err_[i] =
        e1 * k1_[i] + e3 * k3_[i] + e4 * k4_[i] + e5 * k5_[i] + e6 * k6_[i];
  }

  // Error normalization
  ep_ep0 = normaliseError(y_, yTemp_, err_);

  // Step size control
  if (ep_ep0 > 1.0) {
    (*h) *= 0.9 * pow(ep_ep0, -0.25);
    rkStepper(y_, dydx0_, x, h);
  } else {
    (*x) += dx;
    (*h) *= 0.9 * pow(ep_ep0, -0.2);
    for (i = 0; i < N_ode; i++) {
      y_[i] = yTemp_[i];
    }
  }
}

template <int Nmax_ODE_RK>
double RKF45<Nmax_ODE_RK>::normaliseError(double y0[], double y[],
                                          double err[]) {
  double maxErr = 0.0;
  double tol;

  for (int i = 0; i < N_ode; i++) {
    tol = absTol_[i] + relTol_[i] * maxAB(fabs(y0[i]), fabs(y[i]));
    maxErr = maxAB(maxErr, fabs(err[i] / tol));
  }

  return maxErr;
}

template <int Nmax_ODE_RK>
double RKF45<Nmax_ODE_RK>::maxAB(double A, double B) {
  return (A > B) ? A : B;
}

} // namespace numerical
} // namespace utils
} // namespace ziji

#endif // RKF45_H