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

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);
void normalize(std::vector<double>& vec, const double norm);
bool didEscape(int geo, double escape, Vec pos);
Vec calcRandVec(double mach_A, int geo, Vec pos);
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
