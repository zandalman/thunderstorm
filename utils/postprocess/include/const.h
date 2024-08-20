#ifndef CONST_H
#define CONST_H

// includes
#include <string>

constexpr size_t num_mech = 6;
constexpr size_t num_elem = 118;

namespace geo_tag {
constexpr int none = 0;
constexpr int plane = 1;
constexpr int cylinder = 2;
}

namespace stat_tag {
constexpr int ener_loss_mech = 0;
constexpr int num_ion_elem = 1;
constexpr int num_sec_ener = 2;
constexpr int ener_loss_time = 3;
}

namespace bin_tag {
constexpr int ener = 0;
constexpr int escape = 1;
constexpr int ener_sec = 2;
constexpr int time = 3;
}

namespace mech_tag {
constexpr int brem = 0;
constexpr int exc = 1;
constexpr int ion = 2;
constexpr int moller = 3;
constexpr int sync = 4;
constexpr int cher = 5;
}

namespace flags {
constexpr int death = 0;
constexpr int scat = 1;
constexpr int brem = 2;
constexpr int exc = 3;
constexpr int ion = 4;
constexpr int moller = 5;
constexpr int turb = 6;
}

namespace constants {
// length [cm]
constexpr double a0 = 5.29177210e-9; // Bohr radius [cm]
constexpr double AA = 1.e-8;         // Angstrom
constexpr double nm = 1.e-7;         // nanometer
constexpr double um = 1.e-4;         // micrometer
constexpr double mm = 1.e-1;         // millimeter
constexpr double m = 1.e2;           // meter
constexpr double km = 1.e5;          // kilometer
constexpr double R_sol = 6.96e10;    // solar radius
constexpr double AU = 1.496e13;      // astronomical unit
constexpr double ly = 9.463e17;      // lightyear
constexpr double pc = 3.086e18;      // parsec
constexpr double kpc = 3.086e21;     // kiloparsec
constexpr double Mpc = 3.086e24;     // megaparsec

// area [cm^2]
constexpr double sig_T = 6.6524587e-25; // Thompson scattering cross section
constexpr double ba = 1.e-24;           // barn

// mass [g]
constexpr double m_e = 9.1093897e-28;  // electron mass
constexpr double m_p = 1.6726231e-24;  // proton mass
constexpr double m_n = 1.6749286e-24; // neutron mass
constexpr double m_H = 1.660539e-24;   // hydrogen mass
constexpr double amu = 1.6605402e-24;  // atomic mass unit
constexpr double M_sol = 1.9891e+33;   // solar mass

// time [s]
constexpr double min = 60;        // minute
constexpr double hr = 3600;       // hour
constexpr double day = 86400;     // day
constexpr double wk = 604800;     // week
constexpr double yr = 3.1536e7;   // year
constexpr double kyr = 3.1536e10; // kiloyear
constexpr double Myr = 3.1536e13; // megayear
constexpr double Gyr = 3.1536e16; // gigayear

// energy [erg]
constexpr double eV = 1.6021766e-12; // electron volt
constexpr double Ry = 2.1798741e-11; // Rydberg constant
constexpr double keV = 1.6021766e-9; // kiloelectron volt
constexpr double MeV = 1.6021766e-6; // megaelectron volt

// velocity [cm/s]
constexpr double c = 2.9979246e+10; // speed of light [cm/s]

// temperature [K]
constexpr double T_sol = 5.780e3; // solar temperature

// dimensionless
constexpr double alpha = 7.29735256e-3; // fine structure constant
constexpr double Z_sol = 0.0134;        // solar metallicity
constexpr double Y_sol = 0.2485;        // solar helium abundance
constexpr double X_sol = 0.7381;        // solar hydrogen abundance
constexpr double N_A = 6.0221367e23;    // Avagadro's number

// other
constexpr double hbar = 1.0545726e-27; // reduced plank c onstant [erg s]
constexpr double h = 6.6260702e-27;    // plank constant [erg s]
constexpr double k_B = 1.3806490e-16;  // boltzmann constant [erg/K]
constexpr double a_rad =
    7.5657233e-15;                  // radiation density constant [erg/cm^3/K^4]
constexpr double e = 4.8032068e-10; // electron charge [esu]
constexpr double G = 6.67408e-08;   // gravitational constant [cm^3/g/s^2]
constexpr double sig_SB =
    5.67051e-5; // stefan-boltzmann constant [erg/cm^2/K^4/s]
constexpr double L_sol = 3.828e+33; // solar luminosity [erg/s]
} // namespace constants

#endif
