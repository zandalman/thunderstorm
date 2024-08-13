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
double interp(double x0, const Vector1d &x, const Vector1d &y, bool do_llim = false, bool do_ulim = false, double llim = 0., double ulim = 0.);
void calcLamDeb(Vector1d ab, double rho, double temp, double ion_state_avg, double &n_i, double &n_e_free, double &lam_deb);
double calcB0(double n, double temp, double beta, double mach_A);
double calcCosThScat(double xi, double ener, const Vector1d &ener_list, const Vector2d &cos_th_arr_list, const Vector2d &cos_th_dist_list);
double calcEnerLoss(double xi, double ener, const Vector1d &ener_list, const Vector2d &ener_loss_arr_list, const Vector2d &ener_loss_dist_list);
double calcSigBturb(double m_i, double q_i, double gam, double beta, double rho, double Bmag, double Bmag_turb, double q, double Lmax);
double calcPowerSync(double m_i, double q_i, double gam, double beta, double Bmag, double cos_alpha);
double calcPowerCher(double beta, double temp, double n_e_free);
void calcOmxMoller(double gam, double beta, double cos_th_cut, double lam_deb, double &prefac, double &omxmin, double &omxmax, double &omxcut);
double calcPowerMoller(double ener, double gam, double beta, double n_e_free, double lam_deb, double cos_th_cut);
double calcSigMoller(double gam, double beta, double lam_deb, double cos_th_cut);
void calcCosThScatEnerLossMoller(double xi, double ener, double gam, double beta, double lam_deb, double cos_th_cut, double &cos_th, double &ener_loss);
void calcStableParam(double L, double mach_A, double &scale, double &rperp_max);
void calcTurbDiff(double xi1, double xi2, double lam, double scale, double &rperp);
double calcIntermittancyTrans(double xi1, double xi2, double lam_intermittancy, double dis, double cos_alpha);
double calcLamIntermittancy(double m_i, double q_i, double gam, double vmag, double Brms, double L, double mach_A, double alpha);
double calcTurbDiffPar(double xi1, double xi2, double L, double mach_A);

#endif
