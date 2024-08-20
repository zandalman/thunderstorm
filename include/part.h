#ifndef PART_H
#define PART_H

/// @brief A structure to represent a particle.
struct Part {
  int id;                   // The particle ID.
  double m_i;               // The paricle mass [g].
  double q_i;               // The particle charge [esu].
  double ener;              // The particle energy [eV].
  double cos_alpha;         // The particle pitch angle.
  double splus;             // The positive distance travelled along the field line [cm].
  double sminus;            // The negative distance travelled along the field line [cm].
  bool alive;               // Whether the particle is alive.

  Part() = default;
  Part(int id_, double m_i_, double q_i_, double ener_, double cos_alpha_);
  double gam() const;
  double beta() const;
  void scat(double xi, double cos_th);
  void loseEner(double ener_loss);
};

#endif
