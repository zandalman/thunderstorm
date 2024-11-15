// includes
#include <vector>
#include <cmath>
#include <algorithm>
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
 * @brief Determine whether a particle has escaped.
 * 
 * @param geo     The geometry tag.
 * @param escape  The escape distance [cm].
 * @param pos     The particle position [cm].
 * @return Whether the particle escaped or not.    
 */
bool didEscape(int geo, double escape, Vec pos) {
  switch ( geo ) {
    case geo_tag::none:
    return false;
    break;
    case geo_tag::plane:
    return pos.z > escape;
    break;
    case geo_tag::sphere:
    return pos.mag() > escape;
    break;
  }
  return false;
}

Vec calcRandVec(double mach_A) {
  Vec rand_vec = Vec(0.0, 0.0, 0.0);
  double cos_th = 2.0 * xi() - 1.0;
  double sin_th = sqrt(1.0 - cos_th*cos_th);
  double phi = 2.0 * M_PI * xi();
  rand_vec.x = mach_A * sin_th * cos(phi);
  rand_vec.y = mach_A * sin_th * sin(phi);
  rand_vec.z = mach_A * cos_th + 1.0;
  return rand_vec.unit();
}

/**
 * @brief Combine statistics between two datasets.
 * 
 * @param size   The size of the statistics vector.
 * @param nB_int The number of elements in set B.
 * @param meanB  The mean statistics in set B.
 * @param M2B    The M2 statistics in set B.
 * @param M3B    The M3 statistics in set B.
 * @param M4B    The M4 statistics in set B.
 * @param nA_int The number of elements in set A.
 * @param meanA  The mean statistics in set A.
 * @param M2A    The M2 statistics in set A.
 * @param M3A    The M3 statistics in set A.
 * @param M4A    The M4 statistics in set A.
 */
void addStat(
  size_t size, 
  int nB_int, 
  const std::vector<double> &meanB, 
  const std::vector<double> &M2B, 
  const std::vector<double> &M3B,
  const std::vector<double> &M4B,
  int &nA_int, 
  std::vector<double> &meanA, 
  std::vector<double> &M2A,
  std::vector<double> &M3A,
  std::vector<double> &M4A
) {
  double nB = static_cast<double>(nB_int);
  double nA = static_cast<double>(nA_int);
  double nAB = nA + nB;
  double delta, delta_nAB, delta_nAB_sq;
  for ( size_t i = 0; i < size; i++ ) {
    delta = meanB[i] - meanA[i];
    delta_nAB = delta / nAB;
    delta_nAB_sq = delta_nAB*delta_nAB;
    M4A[i] += M4B[i] + delta * delta_nAB * delta_nAB_sq * nA * nB * (nA*nA - nA*nB + nB*nB) \
              + 6.0 * delta_nAB_sq * (nA*nA * M2B[i] + nB*nB * M2A[i]) \
              + 4.0 * delta_nAB * (nA * M3B[i] - nB * M3A[i]);
    M3A[i] += M3B[i] + delta * delta_nAB_sq * nA * nB * (nA - nB) \
              + 3.0 * delta_nAB * (nA * M2B[i] - nB * M2A[i]);
    M2A[i] += M2B[i] + delta * delta_nAB * nA * nB;
    meanA[i] = (nA * meanA[i] + nB * meanB[i]) / nAB;
  }
  nA_int += nB_int;
}

/**
 * @brief Compute central moments from mean, M2, M3, and M4 statistics.
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
  int n_int,
  const std::vector<double> &M2,
  const std::vector<double> &M3,
  const std::vector<double> &M4,
  std::vector<double> &var,
  std::vector<double> &skew,
  std::vector<double> &kurt
) {
  double n = static_cast<double>(n_int);
  double nm1 = n - 1.0;
  for ( size_t i = 0; i < size; i++ ) {
    var[i] = M2[i] / nm1;
    if ( M2[i] == 0 ) {
      skew[i] = 0.0;
      kurt[i] = 0.0;
    } else {
      skew[i] = n * sqrt(nm1) / (n - 2.0) * M3[i] / pow(M2[i], 1.5);
      kurt[i] = n * (n + 1.0) * nm1 * M4[i] / ((n - 2.0) * (n - 3.0) * M2[i]*M2[i]) \
                - 3.0 * nm1*nm1 / ((n - 2.0) * (n - 3.0));
    }
  }
}
