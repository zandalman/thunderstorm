#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <iostream>
#include <vector>
#include <string>

// headers
#include "vec.h"
#include "random.h"
#include "const.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

// table of kernel flags
constexpr int ker_tab[2][2][2] = {
  {
    {ker_tag::loc, ker_tag::perp},
    {ker_tag::perp, ker_tag::cor2}
  },
  {
    {ker_tag::par, ker_tag::cor1},
    {ker_tag::cor1, ker_tag::cor3}
  }
};

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);
void normalize(std::vector<double>& vec, const double norm);
Vec calcRandVec(double mach_A, bool super);
void addStat(
  size_t size, 
  const int nmom,
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
  int nmom,
  int n_int,
  const std::vector<double> &M2,
  const std::vector<double> &M3,
  const std::vector<double> &M4,
  std::vector<double> &var,
  std::vector<double> &skew,
  std::vector<double> &kurt
);
inline int getCellIdx(double x) noexcept;
int calcKer(Vec pos, double dx, int ndim);

#endif
