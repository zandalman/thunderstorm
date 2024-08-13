// includes
#include <cmath>
#include <algorithm>

// headers
#include "const.h"
#include "part.h"
#include "random.h"

/// @brief A constructor for the Part structure.
Part::Part(int id_, double m_i_, double q_i_, double ener_, double cos_alpha_)
  : id(id_)                   // The particle ID.
  , m_i(m_i_)                 // The paricle mass [g].
  , q_i(q_i_)                 // The particle charge [esu].
  , ener(ener_)               // The particle energy [eV].
  , cos_alpha(cos_alpha_)     // The particle pitch angle.
  , splus(0.)                 // The positive distance travelled along the field line [cm].
  , sminus(0.)                // The negative distance travelled along the field line [cm].
  , alive(true)               // Whether the particle is alive.
  {}

/// @brief Compute the particle Lorentz factor.
double Part::gam() const {
  return 1.0 + ener * constants::eV / (m_i * constants::c * constants::c);
}

/// @brief Compute the particle velocity, relative to the speed of light.
double Part::beta() const {
  return sqrt(1.0 - 1.0 / (gam()*gam()));
}

/**
 * @brief Scatter the particle by a given angle.
 * 
 * @param xi     A random number.
 * @param cos_th The cosine of the scattering angle.
*/
void Part::scat(double xi, double cos_th) {
  double sin_alpha = sqrt(1 - cos_alpha*cos_alpha);
  double sin_th = sqrt(1 - cos_th*cos_th);
  double phi = 2.0*M_PI * xi;
  cos_alpha = cos_alpha * cos_th + sin_alpha * sin_th * cos(phi);
}

/**
 * @brief Remove energy from the particle.
 * 
 * @param ener_loss The energy lost by the particle [eV].
*/
void Part::loseEner(double ener_loss) {
  ener = std::max(0., ener - ener_loss);
}
