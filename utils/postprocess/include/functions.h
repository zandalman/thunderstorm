#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <iostream>
#include <vector>

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);
void normalize(std::vector<double>& vec, const double norm);
double arccothMod(double x);
void calcParPerpFac(double mach_A_lpar, double &fac_par, double &fac_perp);
void calcTransport(double mach_A, double L, double s, double dsplus, double dsminus, double &rpar, double &varpar, double &varperp);
void calcFesc(int geo, double escape, double rpar, double varpar, double varperp, double &fesc);
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
  const std::vector1d<double> &M2,
  const std::vector1d<double> &M3,
  const std::vector1d<double> &M4,
  std::vector1d<double> &var,
  std::vector1d<double> &skew,
  std::vector1d<double> &kurt
)

#endif
