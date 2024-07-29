// includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

// headers
#include "functions.h"
#include "part.h"
#include "random.h"

/**
 * @brief Find the index to insert a value into a list.
 * 
 * @param x0     The value.
 * @param x_list The list.
 * @return The index.
*/
int findIdx(double x0, Vector1d x_list) {
  auto it = std::lower_bound(x_list.begin(), x_list.end(), x0);
  return it - x_list.begin();
}

/**
 * @brief Compute the integral of numerical data.
 * 
 * Compute the integral of numerical data using the trapezoidal rule.
 *
 * @param x A vector of x values.
 * @param y A vector of y values.
 * @return The integral of the numerical data.
 */
double integrate(const Vector1d &x, const Vector1d &y) {
  double integral = 0.0;
  for (size_t i = 1; i < x.size(); i++) {
    integral += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  }
  return integral;
}

/**
 * @brief Compute the cumulative distribution function of a variable.
 *
 * @param x       A vector of variable values.
 * @param pdf     A vector of distribution function values.
 * @param initial The initial value of the cumulative distribution function.
 * @return A vector of cumulative distribution functions values.
 */
Vector1d calcCDF(const Vector1d &x, const Vector1d &dist, double initial) {
  Vector1d cdf;
  double integral = integrate(x, dist), running_integral = 0.0;
  cdf.push_back(initial);
  for (size_t i = 1; i < x.size(); i++) {
    running_integral += 0.5 * (dist[i] + dist[i - 1]) * (x[i] - x[i - 1]);
    cdf.push_back(running_integral / integral);
  }
  return cdf;
}

/**
 * @brief Interpolate from numerical data.
 *
 * @param x0      The x value to interpolate.
 * @param x       A vector of x values.
 * @param y       A vector of y values.
 * @param do_llim Use a lower limit on the interpolated value.
 * @param do_ulim Use an upper limit on the interpolated value.
 * @param llim    The lower limit on the interpolated value.
 * @param ulim    The upper limit on the interpolated value.
 * @return The interpolated value.
 */
double interp(double x0, const Vector1d &x, const Vector1d &y, bool do_llim, bool do_ulim, double llim, double ulim) {
  // find left interpolation point
  size_t i = 0;
  if (x0 >= x[x.size() - 2]) {
    i = x.size() - 2;
  } else {
    while (x0 > x[i + 1]) i++;
  }
  // do interpolatation
  double xL = x[i], yL = y[i], xR = x[i + 1], yR = y[i + 1];
  double dydx = (yR - yL) / (xR - xL);
  double y0 = yL + dydx * (x0 - xL);
  // apply limits
  if (do_llim && y0 < llim) {
    y0 = llim;}
  if (do_ulim && y0 > ulim) {
    y0 = ulim;}
  return y0;
}

/**
 * @brief Compute the ion and free electron number densities and the Debye length.
 * 
 * @param ab            The vector of elemental abundances.
 * @param rho           The density [g/cc].
 * @param temp          The temperature [K].
 * @param ion_state_avg The average ionization state.
 * @param n_i           The ion number density [1/cc].
 * @param n_e_free      The free electron number density [1/cc].
 * @param lam_deb       The Debye length [cm].
 */
void calcLamDeb(Vector1d ab, double rho, double temp, double ion_state_avg, double &n_i, double &n_e_free, double &lam_deb) {
  n_i = 0.;
  for ( size_t i = 0; i < ab.size(); i++ ) { n_i += rho * constants::N_A * ab[i]; }
  n_e_free = n_i * ion_state_avg;
  lam_deb  = sqrt(constants::k_B * temp / (4*M_PI * constants::e*constants::e * n_i * ion_state_avg * (ion_state_avg + 1.)));
}

/**
 * @brief Sample the cosine of the scattering angle (theta).
 *
 * @param xi               A random number.
 * @param ener             The energy at which to sample.
 * @param ener_list        A vector of energies [eV].
 * @param cos_th_arr_list  A vector of vectors of cosine theta values at each energy.
 * @param cos_th_dist_list A vector of vectors of cosine theta distribution function values at each energy.
 * @return The sampled cosine theta value.
 */
double calcCosThScat(double xi, double ener, const Vector1d &ener_list, const Vector2d &cos_th_arr_list, const Vector2d &cos_th_dist_list) {
  Vector1d cos_th_arr;
  for (size_t i = 0; i < ener_list.size(); i++) {
    Vector1d cos_th_cdf = calcCDF(cos_th_arr_list[i], cos_th_dist_list[i]);
    cos_th_arr.push_back(
        interp(xi, cos_th_cdf, cos_th_arr_list[i], true, true, -1.0, 1.0));
  }
  return interp(ener, ener_list, cos_th_arr, true, true, -1.0, 1.0);
}

/**
 * @brief Sample the energy loss.
 *
 * @param xi                  A random number.
 * @param ener                The energy at which to sample.
 * @param ener_list           A vector of energies [eV].
 * @param ener_loss_arr_list  A vector of vectors of energy loss values at each energy [eV].
 * @param ener_loss_dist_list A vector of vectors of energy loss distribution function values at each energy [eV].
 * @return The sampled energy loss value.
 */
double calcEnerLoss(double xi, double ener, const Vector1d &ener_list, const Vector2d &ener_loss_arr_list, const Vector2d &ener_loss_dist_list) {
  Vector1d ener_loss_arr;
  for (size_t i = 0; i < ener_list.size(); i++) {
    Vector1d ener_loss_cdf =
        calcCDF(ener_loss_arr_list[i], ener_loss_dist_list[i]);
    ener_loss_arr.push_back(interp(xi, ener_loss_cdf, ener_loss_arr_list[i]));
  }
  return interp(ener, ener_list, ener_loss_arr);
}

/**
 * @brief Compute the synchrotron power.
 * 
 * @param q_i       The particle charge [esu].
 * @param gam       The particle Lorentz factor.
 * @param beta      The particle velocity, relative to the speed of light.
 * @param cos_alpha The cosine of the pitch angle.
 * @param Bmag      The amplitude of the magnetic field [G].
 * @param m_i       The particle mass [g].
 * @return The synchrotron power.
*/
double calcPowerSync(double m_i, double q_i, double gam, double beta, double Bmag, double cos_alpha) {
  return 2./3. * q_i*q_i*q_i*q_i * gam*gam * beta*beta
         * (1 - cos_alpha*cos_alpha) * Bmag*Bmag
         / (m_i*m_i * constants::c*constants::c*constants::c) / constants::eV;
}

/**
 * @brief Compute the Cherenkov power.
 * 
 * When a charged particle moves faster than the plasma phase velocity,
 * it generates perturbations in the electric field (the wakefield)
 * which radiate energy, analogous to Cherenkov radiation.
 * 
 * @param gam      The particle Lorentz factor.
 * @param beta     The particle velocity, relative to the speed of light.
 * @param temp     The background plasma temperature [K].
 * @param n_e_free The free electron density [1/cc].
 * @return The Cherenkov power.
 */
double calcPowerCher(double beta, double temp, double n_e_free) {
  double vel_th = sqrt(3. * constants::k_B * temp / (2. * constants::m_e)); // thermal velocity
  double omega_p = sqrt(4.*M_PI * n_e_free * constants::e*constants::e / constants::m_e); // plasma frequency
  bool subthermal = beta*beta * constants::c*constants::c < 2. * vel_th*vel_th;
  return subthermal ? 0. : constants::e*constants::e * omega_p*omega_p / (2. * beta * constants::c) * log(beta*beta * constants::c*constants::c / (vel_th*vel_th) - 1.) / constants::eV;
}

/**
 * @brief Compute one minus the cosine of the maximum, minimum, and cutoff Moller scattering angle in the CM frame.
 * 
 * @param gam        The particle Lorentz factor.
 * @param beta       The particle velocity, relative to the speed of light.
 * @param cos_th_cut The cosine of the cutoff scattering angle in the lab frame.
 * @param lam_deb    The Debye length [cm].
 * @param prefac     The prefactor used in Moller scattering calculations.
 * @param omxmin     One minus the cosine of the minimum scattering angle in the CM frame.
 * @param omxmax     One minus the cosine of the maximum scattering angle in the CM frame.
 * @param omxcut     One minus the cosine of the cutoff scattering angle in the CM frame.
 */
void calcOmxMoller(double gam, double beta, double cos_th_cut, double lam_deb, double &prefac, double &omxmin, double &omxmax, double &omxcut) {
  prefac = 4*M_PI * constants::e*constants::e*constants::e*constants::e / (constants::m_e*constants::m_e * constants::c*constants::c*constants::c*constants::c * beta*beta) * (gam + 1) / (gam*gam);
  double bmin = constants::h * constants::c / (gam * constants::m_e * beta * constants::c*constants::c);
  double bmax = lam_deb;
  omxmin = prefac * M_PI / (bmin*bmin);
  omxmax = prefac * M_PI / (bmax*bmax);
  omxcut = 2. * (gam + 1.) * (1. - cos_th_cut*cos_th_cut) / (2. + (gam - 1.) * (1. - cos_th_cut*cos_th_cut));
  omxcut = std::min(omxmin, omxcut);
}

/**
 * @brief Compute the small-angle Moller power.
 * 
 * @param gam        The particle Lorentz factor.
 * @param beta       The particle velocity, relative to the speed of light.
 * @param cos_th_cut The cosine of the cutoff scattering angle in the lab frame.
 * @param lam_deb    The Debye length [cm].
 * @param ener       The kinetic energy of the particle [eV].
 * @param n_e_free   The number density of free electrons [1/cc].
 * @return The small-angle Moller power.
 */
double calcPowerMoller(double ener, double gam, double beta, double n_e_free, double lam_deb, double cos_th_cut) {
  double prefac, omxmin, omxmax, omxcut;
  calcOmxMoller(gam, beta, cos_th_cut, lam_deb, prefac, omxmin, omxmax, omxcut);
  return prefac * ener * log(omxcut/omxmax) / 2. * n_e_free * beta * constants::c;
}

/**
 * @brief Compute the cross section of large-angle Moller scattering.
 * 
 * @param gam        The particle Lorentz factor.
 * @param beta       The particle velocity, relative to the speed of light.
 * @param cos_th_cut The cosine of the cutoff scattering angle in the lab frame.
 * @param lam_deb    The Debye length [cm].
 * @return The cross section of large-angle Moller scattering.
 */
double calcSigMoller(double gam, double beta, double lam_deb, double cos_th_cut) {
  double prefac, omxmin, omxmax, omxcut;
  calcOmxMoller(gam, beta, cos_th_cut, lam_deb, prefac, omxmin, omxmax, omxcut);
  return omxcut >= omxmin ? 0. : prefac * (1. / omxcut - 1. / omxmin);
}

/** 
 * @brief Sample the cosine of the scattering angle (theta) and the energy loss for large-angle Moller scattering.
 * 
 * @param xi         A random number.
 * @param gam        The particle Lorentz factor.
 * @param beta       The particle velocity, relative to the speed of light.
 * @param cos_th_cut The cosine of the cutoff scattering angle in the lab frame.
 * @param lam_deb    The Debye length [cm].
 * @param ener       The particle kinetic energy [eV].
 * @param cos_th     The sampled cosine theta value.
 * @param ener_loss  The sampled energy loss value [eV].
 */
void calcCosThScatEnerLossMoller(double xi, double ener, double gam, double beta, double lam_deb, double cos_th_cut, double &cos_th, double &ener_loss) {
  double prefac, omxmin, omxmax, omxcut;
  calcOmxMoller(gam, beta, cos_th_cut, lam_deb, prefac, omxmin, omxmax, omxcut);
  double sig = prefac * (1. / omxcut - 1. / omxmin);
  double omx = 1. / (sig * xi / prefac + 1. / omxmin);
  cos_th = sqrt((2. - omx) * (1. + gam) / (2. * (1 + gam) + omx * (1. - gam)));
  ener_loss = ener * omx / 2.;
}

/**
 * @brief Calculate stable distribution parameters from MHD turbulence parameters.
 * 
 * @param L         The injection scale of the turbulence [cm].
 * @param mach_A    The Alfven Mach number.
 * @param scale     The scale parameter of the stable distribution, normalized by the turbulent mean free path to the power 3/2 [cm^(-1/2)].
 * @param rperp_max The truncation parameter of the stable distribution [cm].
 */
void calcStableParam(double L, double mach_A, double &scale, double &rperp_max) {
  double sig0 = 0.25; // width of stable distribution for scale = 1.0
  double A = tgamma(5.0/3.0) * 0.1875 * sqrt(3.0) / M_PI;
  if ( mach_A < 1.0 ) { // sub-Alfvenic turbulence
    scale = sig0 / sqrt(L) * mach_A*mach_A;
    rperp_max = pow(2.0 * A, -0.75) / sqrt(sig0) * L * mach_A*mach_A*mach_A*mach_A*mach_A;
  } else { // super-Alfvenic turbulence
    scale = sig0 / sqrt(L) * pow(mach_A, 1.5);
    rperp_max = pow(2.0 * A, -0.75) / sqrt(sig0) * L / (mach_A*mach_A*mach_A);
  }
}

/**
 * @brief Calculate the displacement perpendicular to the mean field due to turbulent diffusion.
 * Sample from a truncated stable distribution with stability parameter 2/3.
 *
 * @param xi1, xi2  Two random numbers.
 * @param lam       The turbulent mean free path.
 * @param scale     The scale parameter of the stable distribution, normalized by the turbulent mean free path to the power 3/2 [cm^(-1/2)].
 * @param rperp_max The truncation parameter of the stable distribution [cm].
 * @return The displacement perpendicular to the mean field due to turbulent diffusion [cm].
 */
double calcTurbDiff(double xi1, double xi2, double lam, double scale, double rperp_max) {
  double rperp = 2.0 * rperp_max;
  double U = M_PI * (xi1 - 0.5);
  double W = -log(xi2);
  while ( fabs(rperp) > rperp_max ) {
    rperp = 2.0 * scale * pow(lam, 1.5) / sqrt(W) * sin(U / 3.0) * pow(2.0 * cos(2.0/3.0 * U) - 1.0, 3.0/2.0);
  }
  return rperp;
}
