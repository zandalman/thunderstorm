// includes
#include <cmath>
#include <random>

// headers
#include "const.h"
#include "part.h"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0.0, 1.0);

Vec::Vec(double x_, double y_, double z_): x(x_), y(y_), z(z_) {}
double Vec::mag() const {
    return pow(x*x + y*y + z*z, 0.5);
}
Vec Vec::unit() const {
    double mag_val = mag();
    return Vec(x/mag_val, y/mag_val, z/mag_val);
}

/**
 * @brief Compute a vector in a random direction.
 * 
 * @param mag The magnitude of the random vector.
 * @return The random vector.
*/
Vec randVec(double mag) {
    double phi = 2*M_PI*dis(gen);
    double cos_th = 2*dis(gen)-1;
    double sin_th = pow(1- cos_th*cos_th, 0.5);
    return Vec(mag*sin(phi)*sin_th, mag*cos(phi)*sin_th, mag*cos_th);
}

Part::Part(double m_i_, double q_i_, double ener_): m_i(m_i_), q_i(q_i_), ener(ener_) {
    pos = Vec(0., 0., 0.);
    vel = randVec(beta()*constants::c);
}
double Part::gam() const {
    return 1 + ener*constants::eV/(m_i*constants::c*constants::c);
}
double Part::beta() const {
    double gam_val = gam();
    return pow(1 - 1/(gam_val*gam_val), 0.5);
}
 