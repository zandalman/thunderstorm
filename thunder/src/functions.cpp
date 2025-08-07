// includes
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <mpi.h>

// headers
#include "functions.h"
#include "const.h"
#include "io.h"
#include "vec.h"
#include "random.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

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
  if ( num == 1 ) {
    list.push_back(vmin);
  } else {
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
 * @brief Combine statistics between two datasets.
 * 
 * @param size   The size of the statistics vector.
 * @param nmom.  The number of moments 
 * @param nB_int The number of elements in set B.
 * @param M1B    The mean statistics in set B.
 * @param M2B    The M2 statistics in set B.
 * @param M3B    The M3 statistics in set B.
 * @param M4B    The M4 statistics in set B.
 * @param nA_int The number of elements in set A.
 * @param M1A    The mean statistics in set A.
 * @param M2A    The M2 statistics in set A.
 * @param M3A    The M3 statistics in set A.
 * @param M4A    The M4 statistics in set A.
 */
void addStat(
  size_t size, 
  const int nmom,
  int nB_int, 
  const std::vector<double> &M1B, 
  const std::vector<double> &M2B, 
  const std::vector<double> &M3B,
  const std::vector<double> &M4B,
  int &nA_int, 
  std::vector<double> &M1A, 
  std::vector<double> &M2A,
  std::vector<double> &M3A,
  std::vector<double> &M4A
) {
  double nB = static_cast<double>(nB_int);
  double nA = static_cast<double>(nA_int);
  double nAB = nA + nB;
  double delta, delta_nAB, delta_nAB_sq;
  for ( size_t i = 0; i < size; i++ ) {
    delta = M1B[i] - M2A[i];
    delta_nAB = delta / nAB;
    delta_nAB_sq = delta_nAB*delta_nAB;
    if ( nmom >= 4 ) {
      M4A[i] += M4B[i] + delta * delta_nAB * delta_nAB_sq * nA * nB * (nA*nA - nA*nB + nB*nB) \
              + 6.0 * delta_nAB_sq * (nA*nA * M2B[i] + nB*nB * M2A[i]) \
              + 4.0 * delta_nAB * (nA * M3B[i] - nB * M3A[i]);
    }
    if ( nmom >= 3 ) {
      M3A[i] += M3B[i] + delta * delta_nAB_sq * nA * nB * (nA - nB) \
              + 3.0 * delta_nAB * (nA * M2B[i] - nB * M2A[i]);
    }
    if ( nmom >= 2 ) {
      M2A[i] += M2B[i] + delta * delta_nAB * nA * nB;
    }
    M1A[i] = (nA * M1A[i] + nB * M1B[i]) / nAB;
  }
  nA_int += nB_int;
}

/**
 * @brief Compute central moments from M1, M2, M3, and M4 statistics.
 * 
 * @param size  The size of the statistics vector
 * @param n_int The number of elements.
 * @param M2    The M2 statistics.
 * @param M3    The M3 statistics.
 * @param M4    The M4 statistics.
 * @param var   The variance.
 * @param skew  The skewness.
 * @param kurt  The kurtosis.
 */
void calcMoment(
  size_t size,
  int nmom,
  int n_int,
  const std::vector<double> &M2,
  const std::vector<double> &M3,
  const std::vector<double> &M4,
  std::vector<double> &var,
  std::vector<double> &skew,
  std::vector<double> &kurt
) {
  double n = static_cast<double>(n_int);
  for ( size_t i = 0; i < size; i++ ) {
    var[i] = M2[i] / n;
    if ( M2[i] == 0.0 ) {
      if ( nmom >= 3 ) skew[i] = 0.0;
      if ( nmom >= 4 ) kurt[i] = 0.0;
    } else {
      if ( nmom >= 3 ) skew[i] = sqrt(n) * M3[i] / pow(M2[i], 1.5);
      if ( nmom >= 4 ) kurt[i] = n * M4[i] / (M2[i]*M2[i]) - 3.0;
    }
  }
}

void calcTransportParam(
  double mach_A,
  double dx,
  double lam_turb,
  double &ell_A,
  double &gam_par,
  double &cut_par,
  double &gam_perp,
  double &cut_perp
) {
  double fac = sqrt(2.0); // fudge factor
  double sig0_par = 2.567665550378655; // inter-quartile range for alpha=1/2 Levy stable distribution
  double sig0_perp = 2.2205158963881875; // inter-quartile range for alpha=2/3 Levy stable distribution
  ell_A = mach_A >= 1.0 ? dx / (mach_A*mach_A*mach_A) : dx / (mach_A*mach_A);
  gam_par = sqrt(1.0/180.0) * lam_turb*lam_turb / (sig0_par * ell_A);
  cut_par = 1.0/6.0 * lam_turb*lam_turb / (fac*fac * sig0_par * gam_par);
  gam_perp = sqrt(1.0/9.0 * lam_turb*lam_turb*lam_turb / ell_A) / sig0_perp;
  cut_perp = sqrt(1.0/(6.0*6.0*6.0)) * lam_turb*lam_turb*lam_turb / (fac*fac*fac * sig0_perp*sig0_perp * gam_perp*gam_perp);
}

/**
 * @brief Get the index of the current cell given the normalized 1D position.
 * 
 * @param x The normalized 1D position.
 * @return The index of the current cell.
 */
inline int getCellIdx(double x) noexcept {
  if ( 0.0 <= x && x <= 1.0  ) return 0;
  if ( -1.0 <= x && x <= 2.0 ) return 1;
  return 1;
}

/**
 * @brief Get the kernel flag given the 3D position.
 * 
 * @param pos  The 3D position.
 * @param dx   The cell size.
 * @param ndim The number of dimensions.
 * @return The kernel flag.
 */
int calcKer(Vec pos, double dx, int ndim) {

  const int idx_z = getCellIdx(pos.z / dx);
  const int idx_x = ndim >= 2 ? getCellIdx(pos.x / dx) : -1;
  const int idx_y = ndim >= 3 ? getCellIdx(pos.y / dx) : -1;

  switch (ndim) {
  case 1:
    if ( idx_z == -1 ) {
      return ker_tag::none;
    }
    return ker_tab[idx_z][0][0];
    break;
  case 2:
    if ( idx_z == -1 || idx_x == -1 ) {
      return ker_tag::none;
    }  
    return ker_tab[idx_z][idx_x][0];
    break;
  case 3:
    if ( idx_z == -1 || idx_x == -1 || idx_y == -1 ) {
      return ker_tag::none;
    }
    return ker_tab[idx_z][idx_x][idx_y];
    break;
  default:
    return ker_tag::none;
  }
}

double calcRpar(const double s, const double ell_A) {
  return s <= ell_A 
    ? (1.0 - 1.0/6.0 * s / ell_A) * s 
    : (pow(s / ell_A, 2.0/3.0) - 1.0/6.0) * ell_A;
}

Vec calcTransportStep(
  const double gam_par,
  const double cut_par,
  const double gam_perp,
  const double cut_perp
) {
  Vec step(0.0, 0.0, 0.0);
  double r_perp, cos_th, sgn;
  step.z = gam_par * rvs_stable_1o2(xi(), xi());
  if ( fabs(step.z) > cut_par ) {
    sgn = step.z > 0.0 ? 1.0 : -1.0;
    step.z = sgn * cut_par * rvs_pareto(xi());
  }
  r_perp = gam_perp * rvs_stable_2o3(xi(), xi());
  if ( fabs(r_perp) > cut_perp ) {
    sgn = r_perp > 0.0 ? 1.0 : -1.0;
    r_perp = sgn * cut_perp * rvs_pareto(xi());
  }
  cos_th = 2.0 * (xi() - 0.5);
  step.x = r_perp * cos_th;
  step.y = r_perp * sqrt(1.0 - cos_th*cos_th);
  return step;
}

