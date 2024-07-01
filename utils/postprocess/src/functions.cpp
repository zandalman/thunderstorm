// includes
#include <vector>
#include <cmath>
#include <algorithm>

// headers
#include "functions.h"

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &value_list) {
  value_list.reserve(num);
  if ( log ) {
    vmin = log10(vmin);
    vmax = log10(vmax);
  }
  double step = (vmax - vmin) / (num - 1);
  for ( size_t i = 0; i < num; i++ ) {
    double value = vmin + i * step;
    if ( log ) value = pow(10.0, value);
    value_list.push_back(value);
  }
}

/**
 * @brief Find the index to insert a value into a list.
 * 
 * @param x0     The value.
 * @param x_list The list.
 * @return The index.
*/
size_t findIdx(double x0, std::vector<double> x_list) {
  auto it = std::lower_bound(x_list.begin(), x_list.end(), x0);
  return it - x_list.begin();
}
