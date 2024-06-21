#ifndef PART_H
#define PART_H

// includes
#include "vec.h"

/// @brief A structure to represent a particle.
struct Part {
  int id;      // The particle ID.
  double m_i;  // The paricle mass [g].
  double q_i;  // The particle charge [esu].
  double ener; // The particle energy [eV].
  Vec pos;     // The particle position [cm].
  Vec vel;     // The particle velocity [cm/s].
  Vec Bvec;    // The particle magnetic field [G].

  Part() = default;
  Part(int id_, double m_i_, double q_i_, double ener_);
  double gam() const;
  double beta() const;
  void scat(double cos_th, Vec A);
  void loseEner(double ener_loss);
  void newBvec(double Bmag_co, double Bmag_turb);
};

#endif
