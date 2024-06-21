#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// includes
#include <vector>

// headers
#include "const.h"
#include "part.h"

// types
typedef std::vector<double> Vector1d;
typedef std::vector<std::vector<double>> Vector2d;

int findIdx(double x0, Vector1d x_list);
double integrate(const Vector1d &x, const Vector1d &y);
Vector1d calcCDF(const Vector1d &x, const Vector1d &dist, double initial = 0.0);
double interp(double x0, const Vector1d &x, const Vector1d &y,
              bool do_llim = false, bool do_ulim = false, double llim = 0., double ulim = 0.);
double calcCosThScat(double xi, double ener, const Vector1d &ener_list,
                     const Vector2d &cos_th_arr_list,
                     const Vector2d &cos_th_dist_list);
double calcEnerLoss(double xi, double ener, const Vector1d &ener_list,
                    const Vector2d &ener_loss_arr_list,
                    const Vector2d &ener_loss_dist_list);
double calcSigBturb(const Part &part, double Bmag, double Bmag_turb, double rho,
                    double q, double Lmax);
double calcPowerSync(double q_i, double gam, double beta, double cos_alpha, double Bmag, double m_i);

#endif
