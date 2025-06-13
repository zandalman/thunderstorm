// includes
#include <random>

// headers
#include "random.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

/// @brief Generate a uniformally distributed random number between 0 and 1. 
double xi() { return dis(gen); }

/**
 * @brief Sample a random variable from a standard normal distribution.
 * 
 * Using the Box-Mueller transform.
 * 
 * @return A random variable sampled from a standard normal distribution.
*/
double sampleNormal() {
  double xi1 = xi();
  double xi2 = xi();
  return sqrt(-2.0 * log(xi1)) * cos(2.0 * M_PI * xi2);
}
