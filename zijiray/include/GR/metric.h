#ifndef GR_METRIC_H
#define GR_METRIC_H

#include "objects/astro_objects.h"
#include <cmath>

namespace ziji {
namespace GR {

class Metric {
private:
  const astro_objects::BlackHole &blackHole;

public:
  // Constructor
  Metric(const astro_objects::BlackHole &bh);

  // Method to compute the metric (g_mu_nu)
  void computeMetric(double r, double th, double (*g)[4]) const;

  // Method to compute the first and second derivatives of the metric
  void computeMetricDerivatives(double r, double th, double (*g)[4],
                                double (*dg)[4], double (*dg2)[4]) const;
};

} // namespace GR
} // namespace ziji

#endif // GR_METRIC_H