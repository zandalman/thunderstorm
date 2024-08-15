#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <iostream>
#include <vector>

void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list);
size_t findIdx(double x0, std::vector<double> x_list);
void normalize(std::vector<double>& vec, const double norm);
double arccothMod(double x);
void calcParPerpFac(double mach_A_lpar, double &fac_par, double &fac_perp);
void calcTransport(double mach_A, double L, double s, double dsplus, double dsminus, double &rpar, double &varpar, double &varperp);
void calcFesc(int geo, double escape, double rpar, double varpar, double varperp, double &fesc);

#endif
