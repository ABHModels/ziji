

double FuncRadiEm(double rscr, double *params) {
  double phscr = params[0];
  double rhit = params[1];

  double r_var, th_var;

  ziji::raytrace::RayGen rayGen(black_hole, disk_geom);
  rayGen.FromScreen(observer.distance, observer.inclination, rscr * cos(phscr),
                    rscr * sin(phscr));

  r_var = rayGen.vars[0];
  th_var = rayGen.vars[1];
  params[2] = rayGen.vars[0];
  params[3] = rayGen.impact_par;
  params[4] = rayGen.vars[3];

  double horizon = black_hole.Horizon();
  double tolr = 2.0 * 1.0e-7 * horizon;
  double epsil = 1.e-7;

  // std::cout << "Radi: rem = " << r_var << "\n";

  if (rscr > 0 && r_var > horizon + tolr) {

    return (r_var - rhit) / rhit;

  } else {
    return NAN;
  }
}

void RadiArray(int Nar, double *rArr, double rmax, double rmin) {
  for (int i = 0; i < Nar; i++) {
    rArr[i] = pow(1.0 * i / (Nar - 1), 3) * (rmax - rmin) + rmin;
  }
}
