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
void calcLamDeb(Vector1d ab, double rho, double temp, double q_avg, double qsq_avg, double &n_i, double &n_e_free, double &lam_deb);
double calcB0(double n, double temp, double beta);
double calcCosThScat(double xi, double ener, const Vector1d &ener_list, const Vector2d &cos_th_arr_list, const Vector2d &cos_th_dist_list);
double calcEnerLoss(double xi, double ener, const Vector1d &ener_list, const Vector2d &ener_loss_arr_list, const Vector2d &ener_loss_dist_list);
double calcSigBturb(double m_i, double q_i, double gam, double beta, double rho, double Bmag, double Bmag_turb, double q, double Lmax);
double calcPowerSync(double m_i, double q_i, double gam, double beta, double B0, double cos_alpha);
double calcPowerCerenkov(double beta, double temp, double n_e_free);
void calcOmxMoller(double gam, double beta, double cos_th_cut, double lam_deb, double &prefac, double &omxmin, double &omxmax, double &omxcut);
void calcPowerDiffMoller(
  double ener, 
  double gam, 
  double beta, 
  double cos_alpha,
  double n_e_free, 
  double lam_deb, 
  double cos_th_cut,
  double &edot,
  double &cos_al_dot_sq
);
void calcDiffMott(
  double gam,
  double beta,
  double cos_alpha,
  double n_i,
  double lam_deb,
  double qsq_avg,
  double mmw,
  size_t Zmax,
  const Vector1d &ab,
  double &cos_al_dot_sq
);
double calcSigMoller(double gam, double beta, double lam_deb, double cos_th_cut);
void calcCosThScatEnerLossMoller(
  double xi, 
  double ener, 
  double gam, 
  double beta, 
  double lam_deb, 
  double cos_th_cut, 
  double &cos_th, 
  double &ener_loss
);
void normIonState(Vector1d &ion_state, double &q_avg, double &qsq_avg);

#endif
