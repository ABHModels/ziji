// clang-format off

#include <fstream>
#include <iostream>
#include <istream>
#include <math.h>
#include <string>
#include <vector>

#include <omp.h>

#include <cstddef>

#include <xtensor/xarray.hpp>
#include <xtensor/xnpy.hpp>

#include "return_rad/isotropic.h"

#include "utils/physics.h"
#include "raytrace/raygen.h"

// clang-format on

ziji::astro_objects::BlackHole black_hole(0., 10.);
ziji::astro_objects::AccretionDiskGeometry disk_geom(1.0, 500.0);

using namespace std;
using namespace xt;
using namespace ziji::utils::physics::relativity;

int main() { return 0; }

extern "C" void RemRed(double spin, double rscr, wchar_t *fname, int Npar,
                       int Narry, int Narrx, double *param) {

  wstring ws(fname);
  // your new String
  string fsname(ws.begin(), ws.end());

  vector<vector<int>> cr_arrx(Narry, vector<int>(0));

  vector<double> Thsc(Narrx, 0);
  vector<double> Phsc(Narry, 0);

  // vector<double> Thscnew(2*Narrx, 0);
  // vector<double> Phscnew(2*Narry, 0);

  double thmax = M_PI / 2 - 1.0e-10;
  double phmax = 2 * M_PI;
  double dth = thmax / (Narrx - 1);
  double dph = phmax / (Narry - 1);

  black_hole.init(spin, 10.);

  disk_geom.outer_radius = 500. + 1.0e-5;
  disk_geom.inner_radius = KerrRms(spin);

  double impact_par, gred;
  double remis;

  int i, j, ch; ///, phint;
  // int phint;
  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, ch, remis)

#pragma omp for
  for (i = 0; i < Narry; i++) {
    ziji::raytrace::RayGen rayGen(black_hole, disk_geom);
    ch = 0;
    Phsc[i] = i * dph;
    for (j = 0; j < Narrx; j++) {
      Thsc[j] = j * dth;
      rayGen.FromDisk(rscr, Thsc[j], Phsc[i]);
      remis = rayGen.vars[0];

      if ((remis > disk_geom.inner_radius && remis < disk_geom.outer_radius &&
           fabs(remis - rscr) > 1.0e-5 && fabs(rayGen.vars[4]) > 1.0e-5)) {

        if (ch == 0) {
          if (j == 0)
            cr_arrx[i].push_back(j);
          else {
            cr_arrx[i].push_back(j - 1);
          }

          ch = 1;
        }
        if (ch == 1)
          if (j == Narrx - 1)
            cr_arrx[i].push_back(j);

      } else if (ch == 1) {
        cr_arrx[i].push_back(j);
        ch = 0;
      }
    }
  }

  vector<int> cr_arry;

  //
  ch = 0;

  for (int i = 0; i < cr_arrx.size(); i++) {
    if (cr_arrx[i].size() > 0) {
      if (ch == 0)
        if (i == 0)
          cr_arry.push_back(i);
        else
          cr_arry.push_back(i - 1);
      ch = 1;

      if (ch == 1)
        if (i == Narry - 1)
          cr_arry.push_back(i);
    }

    else if (ch == 1) {
      cr_arry.push_back(i);
      ch = 0;
    }
  }

  int min_th = Narrx;
  for (int i = cr_arry[0]; i < cr_arry[1] + 1; i++) {
    if (cr_arrx[i].size() > 0) {

      if (cr_arrx[i][0] < min_th)
        min_th = cr_arrx[i][0];
    }
  }

  dph = (Phsc[cr_arry[1]] - Phsc[cr_arry[0]]) / (2 * Narry - 1);
  dth = (Thsc[Narrx - 1] - Thsc[min_th]) / (2 * Narrx - 1);

  double dlamb, lamb0, lambmax;

  lamb0 = lambda(Thsc[min_th]);
  lambmax = lambda(Thsc[Narrx - 1]);
  dlamb = (lambmax - lamb0) / (2 * Narrx - 1);

  // cout<< dlamb <<endl;

  xt::xarray<double> redremn = xt::zeros<double>({2 * Narry, 2 * Narrx, 2});

  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, ch, remis)

#pragma omp for
  for (i = 0; i < 2 * Narry; i++) {
    ziji::raytrace::RayGen rayGen(black_hole, disk_geom);

    for (j = 0; j < 2 * Narrx; j++) {

      rayGen.FromDisk(rscr, Thetlam(j * dlamb + lamb0),
                      i * dph + Phsc[cr_arry[0]]);
      remis = rayGen.vars[0];

      if (remis < disk_geom.inner_radius) {
        redremn(i, j, 0) = -1.;
      } else if ((remis < disk_geom.outer_radius &&
                  fabs(remis - rscr) > 1.0e-5 &&
                  fabs(rayGen.vars[4]) > 1.0e-5)) {
        redremn(i, j, 1) =
            RedshiftDisk_Disk(remis, rscr, rayGen.impact_par, black_hole);
        redremn(i, j, 0) = remis;
      }
    }
  }

  param[0] = dph;
  param[1] = dlamb;
  param[2] = Phsc[cr_arry[0]];
  param[3] = Phsc[cr_arry[1]];
  param[4] = lamb0;
  param[5] = lambmax;

  xt ::dump_npy(fsname, redremn);
  // xt :: dump_npy("data/Xsc_Ysc_npy/Xsc.npy", Xscn);
  // xt :: dump_npy("data/Xsc_Ysc_npy/Ysc.npy", Yscn);
}
