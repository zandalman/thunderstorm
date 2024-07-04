#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <iostream>
#include <vector>

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);

#endif
