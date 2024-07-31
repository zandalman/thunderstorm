// includes
#include <cmath>
#include <algorithm>

// headers
#include "const.h"
#include "vec.h"
#include "part.h"
#include "random.h"

/// @brief A constructor for the Part structure.
Part::Part(int id_, double m_i_, double q_i_, double ener_, double B0_)
  : id(id_)                   // The particle ID.
  , m_i(m_i_)                 // The paricle mass [g].
  , q_i(q_i_)                 // The particle charge [esu].
  , ener(ener_)               // The particle energy [eV].
  , pos(Vec(0., 0., 0.))      // The particle position [cm].
  , vel(Vec(0., 0., 0.))      // The particle velocity [cm/s].
  , Bvec(Vec(0., 0., B0_))    // The particle magnetic field [G].
  , alive(true)               // Whether the particle is alive.
  , flag_turb_diff(false)     // Whether to account for turbulent diffusion in transport.
  , flag_intermittancy(false) // Whether to transport using intermittancy.
  , lam_intermittancy(0.0)    // The mean free path to an interaction with a region of high magnetic field curvature.
  {
    vel = randVec(beta() * constants::c);
  }

/// @brief Compute the particle Lorentz factor.
double Part::gam() const {
  return 1.0 + ener * constants::eV / (m_i * constants::c * constants::c);
}

/// @brief Compute the particle velocity, relative to the speed of light.
double Part::beta() const {
  return sqrt(1.0 - 1.0 / (gam()*gam()));
}

/// @brief Pitch angle cosine.
double Part::cos_alpha() const {
  double cos_alpha = dot(vel.unit(), Bvec.unit());
  return cos_alpha;
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
  vel = rotate(vel, A, cos(2.0 * M_PI * xi()));
}

/**
 * @brief Remove energy from the particle.
 * 
 * @param ener_loss The energy lost by the particle [eV].
*/
void Part::loseEner(double ener_loss) {
  ener = std::max(0., ener - ener_loss);
  vel = beta() * constants::c * vel.unit();
}
