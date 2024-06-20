// includes
#include <vector>
#include <cmath>

// headers
#include "part.h"
#include "functions.h"

/**
 * @brief Compute the integral of numerical data.
 * 
 * @param x A vector of x-coordinate values.
 * @param y A vector of y-coordinate values.
 * @return The integral of the numerical data.
*/
double integrate(const Vector1d& x, const Vector1d& y) {
    if ( x.size() != y.size() ) throw std::invalid_argument("x and y must have the same size.");
    double integral = 0.0;
    for (size_t i = 1; i < x.size(); i++) {
        integral += 0.5 * (y[i] + y[i-1]) * (x[i] - x[i-1]);
    }
    return integral;
}

/**
 * @brief Compute the cumulative distribution function of a variable.
 * 
 * @param x A vector of variable values.
 * @param pdf A vector of distribution function values.
 * @param initial The initial value of the cumulative distribution function.
 * @return A vector of cumulative distribution functions values.
*/
Vector1d calcCDF(const Vector1d& x, const Vector1d& dist, double initial) {
    if ( x.size() != dist.size() ) throw std::invalid_argument("x and dist must have the same size.");
    Vector1d cdf;
    double integral = integrate(x, dist), running_integral = 0.0;
    cdf.push_back(initial);
    for (size_t i = 1; i < x.size(); i++) {
        running_integral += 0.5 * (dist[i] + dist[i-1]) * (x[i] - x[i-1]);
        cdf.push_back(running_integral / integral);
    }
    return cdf;
}

/**
 * @brief Interpolate from numerical data.
 * 
 * @param x0 The x-coordinate to interpolate.
 * @param x A vector of x-coordinate values.
 * @param y A vector of y-coordinate values.
 * @param do_llim Use a lower limit on the interpolated value.
 * @param do_ulim Use an upper limit on the interpolated value.
 * @param llim The lower limit on the interpolated value.
 * @param ulim The upper limit on the interpolated value.
 * @return The interpolated value.
*/
double interp(double x0, const Vector1d& x, const Vector1d& y, bool do_llim, bool do_ulim, double llim, double ulim) {
    if ( x.size() != y.size() ) throw std::invalid_argument("x and y must have the same size.");
    // find left interpolation point
    int i = 0;
    if ( x0 >= x[x.size() - 2] ) {
        i = x.size() - 2;
    } else {
        while ( x0 > x[i+1] ) i++;
    }
    // do interpolatation
    double xL = x[i], yL = y[i], xR = x[i+1], yR = y[i+1];
    double dydx = (yR - yL) / (xR - xL);
    double y0 = yL + dydx * (x0 - xL);
    // apply limits
    if ( do_llim && y0 < llim ) y0 = llim;
    if ( do_ulim && y0 > ulim ) y0 = ulim;
    return y0;
}

/**
 * @brief Sample the cosine of the scattering angle (theta).
 * 
 * @param xi A random number.
 * @param ener The energy at which to sample.
 * @param ener_list A vector of energies.
 * @param cos_th_arr_list A vector of vectors of cosine theta values at each energy.
 * @param cos_th_dist_list A vector of vectors of cosine theta distribution function values at each energy.
 * @return The sampled cosine theta value.
*/
double calcCosThScat(double xi, double ener, const Vector1d& ener_list, const Vector2d& cos_th_arr_list, const Vector2d& cos_th_dist_list) {
    if ( ener_list.size() != cos_th_arr_list.size() ) throw std::invalid_argument("ener_list and cos_th_arr_list must have the same size.");
    if ( ener_list.size() != cos_th_dist_list.size() ) throw std::invalid_argument("ener_list and cos_th_dist_list must have the same size.");
    Vector1d cos_th_arr;
    for (size_t i=0; i<ener_list.size(); i++) {
        Vector1d cos_th_cdf = calcCDF(cos_th_arr_list[i], cos_th_dist_list[i]);
        cos_th_arr.push_back(interp(xi, cos_th_cdf, cos_th_arr_list[i], true, true, -1.0, 1.0));
    }
    return interp(ener, ener_list, cos_th_arr, true, true, -1.0, 1.0);
}

/**
 * @brief Sample the energy loss.
 * 
 * @param xi A random number.
 * @param ener The energy at which to sample.
 * @param ener_list A vector of energies.
 * @param ener_loss_arr_list A vector of vectors of energy loss values at each energy.
 * @param ener_loss_dist_list A vector of vectors of energy loss distribution function values at each energy.
 * @return The sampled energy loss value.
*/
double calcEnerLoss(double xi, double ener, const Vector1d& ener_list, const Vector2d& ener_loss_arr_list, const Vector2d& ener_loss_dist_list) {
    if ( ener_list.size() != ener_loss_arr_list.size() ) throw std::invalid_argument("ener_list and ener_loss_arr_list must have the same size.");
    if ( ener_list.size() != ener_loss_dist_list.size() ) throw std::invalid_argument("ener_list and ener_loss_dist_list must have the same size.");
    Vector1d ener_loss_arr;
    for (size_t i=0; i<ener_list.size(); i++) {
        Vector1d ener_loss_cdf = calcCDF(ener_loss_arr_list[i], ener_loss_dist_list[i]);
        ener_loss_arr.push_back(interp(xi, ener_loss_cdf, ener_loss_arr_list[i]));
    }
    return interp(ener, ener_list, ener_loss_arr);
}

/**
 * @brief Compute the effective cross section for interactons between a particle and a turbulent magnetic field.
 * 
 * @param part The particle struct.
 * @param Bmag The amplitude of the total magnetic field.
 * @param Bmag_turb The amplitude of the turbulent magneitc field.
 * @param rho The density of the background plasma.
 * @param q The power law exponent of the magnetic turbulence spectrum.
 * @param Lmax The largest scale of magnetic turbulence [cm].
 * @return The effective cross section [cm^2].
*/
double calcSigBturb(const Part& part, double Bmag, double Bmag_turb, double rho, double q, double Lmax) {
    double beta_val = part.beta(), gam_val = part.gam();
    double beta_A = Bmag / pow(4*M_PI*rho*constants::c*constants::c, 0.5); // Alfven speed relative to the speed of light
    double func_beta_A = (1 - pow(beta_A/beta_val, 2-q)) / (2-q) - (1 - pow(beta_A/part.beta(), 4-q)) / (4-q);
    double Om = part.q_i*Bmag/(part.m_i*constants::c); // particle gyro-frequency
    double kmin = constants::c/(Om*Lmax); // minimum wavenumber of magnetic turbulence spectrum
    double fturb = pow(Bmag_turb/Bmag, 2.0); // fraction of magnetic energy in the turbulent field
    double lam_Bturb = constants::c*pow(beta_val*gam_val, 2-q)*func_beta_A * 2.0 / (M_PI*(q-1.0)*fturb*constants::c*kmin) * pow(constants::c*kmin/Om, 2.0-q);
    return 1/(lam_Bturb*rho*constants::N_A);
}

Vec rotateVec(Vec v1, Vec v2) {
    Vec v3(0., 0., 0.);
    return v3;
}
