#ifndef VEC_H
#define VEC_H

/// @brief A structure to represent a 3D vector.
struct Vec {
  double x, y, z; // The vector components.

  Vec() = default;
  Vec(double x_, double y_, double z_);
  double mag() const;
  Vec unit() const;
  Vec operator +(const Vec& other) const;
  Vec operator *(double k) const;
  friend Vec operator *(double k, const Vec& self);
};

Vec cross(Vec A, Vec B);
double dot(Vec A, Vec B);
Vec rotate(Vec A, Vec B, double cos_th);
Vec randVec(double mag = 1.0);

#endif
