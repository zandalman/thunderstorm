#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <iostream>
#include <vector>

// headers
#include "vec.h"
#include "random.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief Miscellaneous parameters
struct MiscParam {
  double rho_sim;  // The density in the simulation [g/cc]
  double ener_min; // The minimum energy in the simulation [eV]
  double inner;    // The thermalization barrier position [scale]
  double outer;    // The escape barrier position [scale]
  double turb;     // The turbulence injection scale [scale]

  MiscParam() = default;
  MiscParam(
    double rho_sim_,
    double ener_min_,
    double inner_,
    double outer_,
    double turb_
  );
};

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);
void normalize(std::vector<double>& vec, const double norm);
Vec calcRandVec(double mach_A);
void addStat(
  size_t size, 
  int nB_int, 
  const std::vector<double> &meanB, 
  const std::vector<double> &M2B, 
  const std::vector<double> &M3B,
  const std::vector<double> &M4B,
  int &nA_int, 
  std::vector<double> &meanA, 
  std::vector<double> &M2A,
  std::vector<double> &M3A,
  std::vector<double> &M4A
);
void calcMoment(
  size_t size,
  int n_int,
  const std::vector<double> &M2,
  const std::vector<double> &M3,
  const std::vector<double> &M4,
  std::vector<double> &var,
  std::vector<double> &skew,
  std::vector<double> &kurt
);

#endif
