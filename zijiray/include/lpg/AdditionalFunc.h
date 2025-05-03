
double FuncRadiEm(double del, double *params) {
  double hlp = params[0];
  double r_var, th_var;
  double rhit = params[1];

  ziji::raytrace::RayGen rayGen(black_hole, disk_geom);
  rayGen.FromLampPost(hlp, del);
  // radi_emis(hlp, del, &impact_par, vars)

  r_var = rayGen.vars[0];
  th_var = rayGen.vars[1];
  params[2] = rayGen.vars[0];
  params[3] = rayGen.impact_par;

  // Bh change

  if (r_var < rhor + 1.0e-4)
    return (r_var * fabs(sin(th_var)) - rhit) / rhit;

  return (r_var - rhit) / rhit; // dr/rhit = (r-rhit)/rhit
}

void RadiArray(int Nar, double *rArr, double rmax, double rmin) {
  for (int i = 0; i < Nar; i++) {
    rArr[i] = pow(1.0 * i / (Nar - 1), 3) * (rmax - rmin) + rmin;
  }
}

void diffArray(int Nar, double *Arr, double *dArr) {

  for (int i = 0; i < Nar - 1; i++) {
    dArr[i] = Arr[i + 1] - Arr[i];
  }

  dArr[Nar - 1] = dArr[Nar - 2];
}

void Save_XY(const char *filename, double *xarray, double *yarray, int size) {

  FILE *out;

  out = fopen(filename, "w");

  // Write the data to the file
  for (int i = 0; i < size; i++) {
    fprintf(out, "%e %e\n", xarray[i], yarray[i]);
  }

  // Close the file
  fclose(out);
}
