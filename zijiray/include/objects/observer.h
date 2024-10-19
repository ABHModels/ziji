#ifndef OBSERVER_H
#define OBSERVER_H

namespace ziji {

// Observer struct
struct Observer {
  double distance;    // in rg
  double inclination; // in degrees
};

// ScreenNumberPhotons struct with constructor
struct ScreenNumberPhotons {
  int Npsc;
  int Nrsc;

  // Constructor that allows initializing Npsc and Nrsc
  ScreenNumberPhotons(int npsc, int nrsc) : Npsc(npsc), Nrsc(nrsc) {}
};

// ReturnNumberPhotons struct with constructor
struct ReturnNumberPhotons {
  int Nth;
  int Nph;

  // Constructor that allows initializing Nth and Nph
  ReturnNumberPhotons(int nth, int nph) : Nth(nth), Nph(nph) {}
};

} // namespace ziji

#endif // OBSERVER_H