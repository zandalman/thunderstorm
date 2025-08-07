// includes
#include <random>
#include <cmath>

// headers
#include "random.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

/**
 * @brief Generate a uniformally distributed random number between 0 and 1.
 */
double xi() { return dis(gen); }

// valgrind friendly version for debugging
// double xi() {
//   static uint32_t seed = 123456789;  // You can set this to any number
//   seed = 1664525 * seed + 1013904223;  // LCG parameters from Numerical Recipes
//   return static_cast<double>(seed) / static_cast<double>(UINT32_MAX);
// }

/**
 * @brief Sample an exponential distribution.
 */
double rvs_exp(double xi) {
  return -log(xi);
}

/**
 * @brief Sample a Pareto distribution with alpha=1.
 */
double rvs_pareto(double xi) {
  return 1.0 / (1.0 - xi);
}

/**
 * @brief Sample a Levy stable distribution with alpha=1/2.
 */
double rvs_stable_1o2(double xi1, double xi2) {
  double U = M_PI * (xi1 - 0.5);
  double W = rvs_exp(xi2);
  return 0.5 / W * tan(U) / cos(U);
}

/**
 * @brief Sample a Levy stable distribution with alpha=1/2.
 */
double rvs_stable_2o3(double xi1, double xi2) {
  double U = M_PI * (xi1 - 0.5);
  double W = rvs_exp(xi2);
  double temp = sqrt(2.0 * cos(2.0/3.0 * U) - 1.0);
  return 2.0 / sqrt(W) * sin(U/3.0) / (temp*temp*temp);
}
