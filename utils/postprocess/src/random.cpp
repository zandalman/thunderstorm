// includes
#include <random>

// headers
#include "random.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

/// @brief Generate a uniformally distributed random number between 0 and 1. 
double xi() { return dis(gen); }
