// includes
#include <cmath>

// headers
#include "const.h"
#include "random.h"
#include "vec.h"

/// @brief A constructor for the Vec structure.
Vec::Vec(double x_, double y_, double z_) 
  : x(x_), y(y_), z(z_) // The vector components.
  {}

/// @brief Compute the magnitude of the vector.
double Vec::mag() const { 
  return pow(x * x + y * y + z * z, 0.5);
}

/// @brief Compute the unit vector.
Vec Vec::unit() const {
  double mag_val = mag();
  return Vec(x / mag_val, y / mag_val, z / mag_val);
}

/// @brief Overload the + operator for vector addition.
Vec Vec::operator +(const Vec& other) const {
  return Vec(x + other.x, y + other.y, z + other.z);
}

/// @brief Overload the * operator for scalar-vector multiplication.
Vec Vec::operator *(double k) const {
  return Vec(x * k, y * k, z * k);
}
Vec operator *(double k, const Vec& self) {
  return self * k;
}

/**
 * @brief Take the dot product of two vectors.
 * 
 * @param A The first vector.
 * @param B The second vector.
 * @return The dot product of the vectors.
*/
double dot(Vec A, Vec B) {
  return A.x * B.x + A.y * B.y + A.z * B.z;
}

/**
 * @brief Take the cross product of two vectors.
 * 
 * @param A The first vector.
 * @param B The second vector.
 * @return The cross product of the vectors.
*/
Vec cross(Vec A, Vec B) {
  return Vec(
    A.y * B.z - A.z * B.y,
    A.z * B.x - A.x * B.z,
    A.x * B.y - A.y * B.x
  );
}

/**
 * @brief Rotate one vector about another vector.
 * 
 * Rotate one vector about another vector by a given angle
 * using the Rodrigues rotation formula.
 * 
 * @param A The vector to rotate.
 * @param B The vector about which to rotate.
 * @param th The cosine of the rotation angle.
 * @return The rotated vector.
*/
Vec rotate(Vec A, Vec B, double cos_th) {
  double sin_th = pow(1 - cos_th*cos_th, 0.5);
  Vec Bhat = B.unit();
  return A * cos_th + cross(Bhat, A) * sin_th + Bhat * dot(Bhat, A) * (1 - cos_th);
}

/**
 * @brief Compute a vector in a random direction.
 *
 * @param mag The magnitude of the random vector.
 * @return The random vector.
 */
Vec randVec(double mag) {
  double phi = 2 * M_PI * xi();
  double cos_th = 2 * xi() - 1;
  double sin_th = pow(1 - cos_th * cos_th, 0.5);
  return Vec(mag * sin(phi) * sin_th, mag * cos(phi) * sin_th, mag * cos_th);
}
