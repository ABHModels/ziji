// clang-format off
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <string>

#include <omp.h>

#include <cstddef>

#include <xtensor/xarray.hpp>
#include <xtensor/xnpy.hpp>



using namespace std;
using namespace xt;


#include "utils/physics.h"
#include "raytrace/raygen.h"
#include "objects/observer.h"

ziji::astro_objects::BlackHole black_hole(0., 10.);
ziji::astro_objects::AccretionDiskGeometry disk_geom(1.0, 500.0);

ziji::Observer observer;

#include "utils/f_root.h"
using namespace ziji::utils::physics::relativity;
using namespace ziji::utils::numerical;

#include "transit/AdditionalFunc.h"

// clang-format on

int main() { return 0; }

void initGlobals(double spin, double incang) {
  observer.inclination = M_PI * incang / 180.;
  observer.distance = 10000.;

  black_hole.init(spin, 10.);

  disk_geom.outer_radius = 500.;
  disk_geom.inner_radius = KerrRms(spin);
}

extern "C" void RadiRedshArray(double spin, double incang, wchar_t *fname,
                               int Npar) {
  initGlobals(spin, incang);

  auto Xsc = load_npy<double>("data/data_cache/transit_ray/Xsc.npy");
  auto Ysc = load_npy<double>("data/data_cache/transit_ray/Ysc.npy");

  wstring ws(fname);
  // your new String
  string fsname(ws.begin(), ws.end());

  auto Xshape = adapt(Xsc.shape());

  int Narry = Xshape(0);
  int Narrx = Xshape(1);
  xt::xarray<double> RadRed = xt::zeros<double>({Narry, Narrx, 2});

  double remis;

  // int Nth = 8;
  int i, j;

  omp_set_num_threads(Npar);
#pragma omp parallel private(i, j, remis)

#pragma omp for
  for (i = 0; i < Narry; i++) {
    ziji::raytrace::RayGen rayGen(black_hole, disk_geom);
    for (j = 0; j < Narrx; j++) {
      // radi_emis(Xsc(i, j), Ysc(i, j), &impact_par, vars);
      rayGen.FromScreen(observer.distance, observer.inclination, Xsc(i, j),
                        Ysc(i, j));

      remis = rayGen.vars[0];

      if ((remis > disk_geom.inner_radius &&
           remis < disk_geom.outer_radius + 1.0e-4)) {
        RadRed(i, j, 0) = remis;
        RadRed(i, j, 1) =
            RedshiftDisk_Screen(remis, rayGen.impact_par, black_hole);

        // printf("%f %f %d %d %f %f \n",Flux(i,j),remis, i,j, Xsc(i,j),
        // Ysc(i,j));
      }
    }
  }

  xt ::dump_npy(fsname, RadRed);
}

extern "C" void rdiskLineScr(int nr_zone, double spin, double incang,
                             wchar_t *fname) {
  initGlobals(spin, incang);

  auto phsc = load_npy<double>("data/data_cache/transit_ray/phsc.npy");

  wstring ws(fname);
  // your new String
  string fsname(ws.begin(), ws.end());

  auto phscshape = adapt(phsc.shape());

  int Narr_ph = phscshape(0);

  xt::xarray<double> Rscr = xt::zeros<double>({Narr_ph});

  double err_radi = 1.0e-7;

  double remis;
  double rscr_Sol;

  double rdisk[2];

  remis = rdisk[nr_zone];

  if (nr_zone == 0)
    remis = disk_geom.inner_radius;
  else
    remis = disk_geom.outer_radius;

  double rem0 = remis * cos(observer.inclination);
  double rem1 = rem0 * 1.1;

  if (remis < 5. * black_hole.Horizon()) {
    rem0 = 2 * remis * cos(observer.inclination);
    rem1 = rem0 * 2;
  }

  double params[4];

  int i;

  for (i = 0; i < Narr_ph; i++) {

    params[0] = phsc(i);
    params[1] = remis;

    rscr_Sol = scMethod(FuncRadiEm, params, rem0, rem1, err_radi, 400);

    if (rscr_Sol < 0) {
      printf("No %d:  rscr = %f;  phsc = %f;  rem=  %f; \n", i, rscr_Sol,
             params[0], remis);
    }
    rem0 = rscr_Sol;
    rem1 = rscr_Sol + rscr_Sol * err_radi * 100.;
    Rscr(i) = rscr_Sol;
    // Rscr_Redsh(i, 1) = redshift(remis, impact_par);
  }

  xt ::dump_npy(fsname, Rscr);
}
