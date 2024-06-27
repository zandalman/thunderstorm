// includes
#include <cmath>
#include <algorithm>

// headers
#include "const.h"
#include "vec.h"
#include "part.h"
#include "random.h"

/// @brief A constructor for the Part structure.
Part::Part(int id_, double m_i_, double q_i_, double ener_)
  : id(id_)               // The particle ID.
  , m_i(m_i_)             // The paricle mass [g].
  , q_i(q_i_)             // The particle charge [esu].
  , ener(ener_)           // The particle energy [eV].
  , pos(Vec(0., 0., 0.))  // The particle position [cm].
  , vel(Vec(0., 0., 0.))  // The particle velocity [cm/s].
  , Bvec(Vec(0., 0., 0.)) // The particle magnetic field [G].
  {
    vel = randVec(beta() * constants::c);
  }

/// @brief Compute the particle Lorentz factor.
double Part::gam() const {
  return 1.0 + ener * constants::eV / (m_i * constants::c * constants::c);
}

/// @brief Compute the particle velocity, relative to the speed of light.
double Part::beta() const {
  double gam_val = gam();
  return sqrt(1.0 - 1.0 / (gam_val * gam_val));
}

/**
 * @brief Scatter the particle by a given angle off a given axis.
 * 
 * Scatter the particle by a given angle off a given axis:
 *  1. Compute a unit vector perpendicular to the scattering axis.
 *  2. Rotate the velocity about the perpendicular vector by the scattering angle.
 *  3. Rotate the velocity about the scattering direction by a random angle.
 * 
 * @param cos_th The cosine of the scattering angle.
 * @param A The scattering axis.
*/
void Part::scat(double cos_th, Vec A) {
  Vec vperp = cross(A, randVec());
  vel = rotate(vel, vperp, cos_th);
  vel = rotate(vel, A, cos(2*M_PI*xi()));

}

/**
 * @brief Remove energy from the particle.
 * 
 * @param ener_loss The energy lost by the particle [eV].
*/
void Part::loseEner(double ener_loss) {
  ener = std::max(0., ener - ener_loss);
}

/**
 * @brief Resample the particle magnetic field.
 * 
 * @param Bmag_co   The amplitude of the coherent magnetic field [G].
 * @param Bmag_turb The amplitude of the turbulent magnetic field [G].
*/
void Part::newBvec(double Bmag_co, double Bmag_turb) {
  Bvec = Bmag_co * Vec(0., 0., 1.) + randVec(Bmag_turb);
}
