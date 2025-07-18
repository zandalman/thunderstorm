// includes
#include <random>

// headers
#include "random.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

// @brief Generate a uniformally distributed random number between 0 and 1. 
// double xi() { return dis(gen); }

// valgrind friendly version
double xi() {
  static uint32_t seed = 123456789;  // You can set this to any number
  seed = 1664525 * seed + 1013904223;  // LCG parameters from Numerical Recipes
  return static_cast<double>(seed) / static_cast<double>(UINT32_MAX);
}
