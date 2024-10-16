// clang-format off
#include <istream>
#include <fstream>
#include <cmath>

#include <string>       // std::string
#include <iostream>     // std::cout
#include <stdio.h>
//#include <omp.h>

#include <cstddef>

#include "utils/physics.h"
#include "lpg/global_var.h"
#include "raytrace/raygen.h"

ziji::astro_objects::BlackHole black_hole(0.,10.);
ziji::astro_objects::AccretionDiskGeometry disk_geom(
      1.0,
      500.0);


#include "utils/f_root.h"

#include "lpg/AdditionalFunc.h"


using namespace ziji::utils::physics::relativity;
using namespace ziji::utils::numerical;

using namespace std;

// clang-format on

// void LPGeom(double spin, double hlp);

int main(int argc, char *argv[]) {
  // double a, height;

  // a = atof(argv[1]);      /* spin parameter */
  // height = atof(argv[2]); /* height of the source */

  // LPGeom(a, height);

  return 0;
}

extern "C" void LPGeom(double spin, double hlp) //, wchar_t* fname)
{

  absTol_r = 1.0e-4;

  double err_radi = 1.0e-5;

  // rhor = 1.0 + sqrt(1.0 - chi * chi);
  black_hole.init(spin, 10.);
  rhor = black_hole.Horizon();
  isco = KerrRms(spin);
  // cout << isco << "   " << rhor << "\n";
  disk_geom.inner_radius = isco;
  disk_geom.outer_radius = 1000;
  // disk_radi_max = 1500.;

  double impact_par;
  double vars[5];
  double remis;

  double radArr[Nradi_prof];
  double radArrHit[Nradi_prof];
  double diffradArr[Nradi_prof];
  double IntenArr[Nradi_prof]; // Inten = sin(del) ddel/dr 1/(Lorens*dA/dr) =
                               // dOmega_obs/dA_obs 1/(2pi)
  double Lor_dA_dr;
  double rmin = isco * (1.0 + err_radi);
  double rmax = 1000.;

  RadiArray(Nradi_prof, radArr, rmax, rmin);
  diffArray(Nradi_prof, radArr, diffradArr);

  double outRoot[3];
  double params[4];

  double del0 = 1.0e-6;
  double del1 = 1.5 * 1.0e-6;
  double delhit, rhit;
  double ddel_dr;

  params[0] = hlp;

  for (int i = 0; i < Nradi_prof; i++) {
    params[1] = radArr[i];
    scMethod2(FuncRadiEm, params, outRoot, del0, del1, err_radi, 100);
    delhit = outRoot[0];
    rhit = (outRoot[1] + 1.0) * radArr[i];
    ddel_dr = 1 / (outRoot[2] * radArr[i]);

    del0 = delhit + ddel_dr * diffradArr[i];
    del1 = delhit + ddel_dr * diffradArr[i] * (1.0 + 1.0e-4);

    Lor_dA_dr = Loren_dA_dr(rhit, isco, black_hole);

    radArrHit[i] = rhit;
    IntenArr[i] = sin(delhit) * ddel_dr / Lor_dA_dr;

    // cout<<radArr[i]<<"   "<<radArrHit[i]<<"  "<<IntenArr[i]<< endl;
  }

  char filename[200];
  snprintf(filename, sizeof(filename),
           "data/data_local/lampost/lp_%.3f_%.3f.dat", spin, hlp);

  Save_XY(filename, radArrHit, IntenArr, Nradi_prof);

  // printf("File: %s  \n", filename);
}
