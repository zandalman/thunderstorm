// includes
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>

// headers
#include "functions.h"
#include "const.h"

/**
 * @brief Make a linearly or logarithmically spaced list of values.
 * 
 * @param vmin The minimum value.
 * @param vmax The maximum value.
 * @param num  The number of values.
 * @param log  Whether the list is logarithmically spaced.
 * @param list The vector to store the list.
 */
void linspace(double vmin, double vmax, size_t num, bool log, std::vector<double> &list) {
  list.reserve(num);
  if ( log ) {
    vmin = log10(vmin);
    vmax = log10(vmax);
  }
  double step = (vmax - vmin) / (num - 1);
  for ( size_t i = 0; i < num; i++ ) {
    double value = vmin + i * step;
    if ( log ) value = pow(10.0, value);
    list.push_back(value);
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

/**
 * @brief Normalize a vector by a constant.
 * 
 * @param vec The vector to normalize.
 * @param norm The constant by which to normalize the vector.
 */
void normalize(std::vector<double>& vec, const double norm) {
  for ( size_t i = 0; i < vec.size(); i++ ) {
    vec[i] /= norm;
  }
}

/**
 * @brief Modified inverse hyperbolic cotangent function.
 * The domain of the inverse hyperbolic cotangent is extended by taking the maximum of x and 1/x.
 * 
 * @param x The argument of the function.
 * @return The modified inverse hyperbolic cotangent.
 */
double arccothMod(double x) {
  double x_mod = fmax(x, 1.0 / x);
  return 0.5 * log((x_mod + 1.0) / (x_mod - 1.0));
}

/**
 * @brief Calculate parallel and perpendicular diffusion factors.
 * 
 * @param mach_A_lpar The parallel scale-dependent Alfven Mach number.
 * @param fac_par     The parallel diffusion factor.
 * @param fac_perp    The perpendicular diffusion factor.
 */
void calcParPerpFac(double mach_A_lpar, double &fac_par, double &fac_perp) {
  if ( mach_A_lpar == 0.0 ) {
    fac_par = 0.; fac_perp = 0.;
  } else if ( mach_A_lpar == 1.0 ) {
    fac_par = 0.5; fac_perp = 0.5;
  } else {
    double mach_A_lpar_sqm1 = mach_A_lpar*mach_A_lpar - 1.0;
    double fac = mach_A_lpar_sqm1*mach_A_lpar_sqm1 * arccothMod(mach_A_lpar) / (4.0 * mach_A_lpar);
    fac_par = fac - (mach_A_lpar_sqm1 - 2.0) / 4.0;
    fac_perp = -fac + (mach_A_lpar_sqm1 + 2.0) / 4.0;
  }
}

/**
 * @brief Calculate the parallel and perpendicular transport.
 * 
 * @param mach_A  The Alfven mach number.
 * @param L       The turbulence injection scale [cm].
 * @param s       The arclength along the field line [cm].
 * @param dsplus  The differential postive arclength along the field line [cm].
 * @param dsminus The differential negative arclength along the field line [cm].
 * @param rpar    The mean parallel displacement [cm].
 * @param sigpar  The standard deviation of the parallel displacement [cm].
 * @param sigperp The standard deviation of the perpendicular displacement [cm].
 */
void calcTransport(double mach_A, double L, double s, double dsplus, double dsminus, double &rpar, double &varpar, double &varperp) {
  double ds, sgn, lcorr;
  double mach_A_lpar, drpar_ds, fac_par, fac_perp;
  // compute ds magnitude and sign
  ds = dsplus + dsminus;
  sgn = dsplus > 0.0 ? 1.0 : -1.0;
  // compute the scale-dependent Mach number
  if ( s >= L ) {
    mach_A_lpar = mach_A;
  } else if ( mach_A < 1 ) {
    mach_A_lpar = mach_A*mach_A * sqrt(s / L);
  } else if ( s >= L / (mach_A*mach_A*mach_A) ) {
    mach_A_lpar = mach_A * pow(s / L, 1.0/3.0);
  } else {
    mach_A_lpar = pow(mach_A, 3.0/2.0) * sqrt(s / L);
  }
  // compute useful quantities
  lcorr = s >= L ? L : s;
  drpar_ds = mach_A_lpar >= 1 ? 2.0/3.0 / mach_A_lpar : 1.0 - mach_A_lpar*mach_A_lpar / 3.0;
  calcParPerpFac(mach_A_lpar, fac_par, fac_perp);
  // increment rpar, sigpar, and sigperp
  rpar += sgn * ds * drpar_ds;
  varpar += fmax(0.0, ds * lcorr * (fac_par - drpar_ds*drpar_ds));
  varperp += fmax(0.0, ds * lcorr * fac_perp);
}

/**
 * @brief Compute the survival fraction.
 * 
 * @param geo     The geometry tag.
 * @param escape  The escape distance [cm].
 * @param rpar    The mean parallel displacement [cm].
 * @param sigpar  The standard deviation of the parallel displacement [cm].
 * @param sigperp The standard deviation of the perpendicular displacement [cm].
 * @param fesc    The survival fraction.
 */
void calcFesc(int geo, double escape, double rpar, double varpar, double varperp, double &fesc) {
  double sigpar = sqrt(varpar);
  double sigperp = sqrt(varperp);
  switch ( geo ) {
    case geo_tag::none:
    fesc = 1.;
    break;
    case geo_tag::plane:
    fesc = 0.5 * (1.0 + erf((escape - rpar) / sqrt(2.0) / sigpar));
    break;
    case geo_tag::cylinder:
    fesc = 0.5 * (1.0 - exp(-escape*escape / (2.0 * sigperp*sigperp))) * (erf((escape - 2.0 * rpar) / (sqrt(8.0) * sigpar)) + erf((escape + 2.0 * rpar) / (sqrt(8.0) * sigpar)));
    break;
  }
}
