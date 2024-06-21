''' Units and physical constants. '''
import numpy as np

# length units
AA  = 1e-8      # Angstroms [cm]
um  = 1e-4      # micron [cm]
m   = 100       # meter [cm]
km  = 1000*m    # kilometer [cm]
AU  = 1.496e13  # astronomical unit [cm]
ly  = 9.463e17  # lightyear [cm]
pc  = 3.086e+18 # parsec [cm]
kpc = 1000*pc   # kiloparsec [cm]
Mpc = 1e6*pc    # megaparsec [cm]

# area units
ba  = 1e-24     # barn [cm^2]

# time units
hr  = 3600    # hour [s]
day = 24*hr   # day [s]
wk  = 7*day   # week [s]
yr  = 365*day # year [s]
kyr = 1e3*yr  # kiloyear [s]
Myr = 1e6*yr  # megayear [s]
Gyr = 1e9*yr  # gigayear [s]

# physical constants
alpha  = 0.00729735256 # fine structure constant
a0     = 5.291772e-9   # Bohr radius [cm]
c      = 2.9979246e+10 # speed of light [cm/s]
h      = 6.6260702e-27 # plank constant [erg s]
hbar   = h/(2*np.pi)   # reduced plank c onstant [erg s]
G      = 6.67408e-08   # gravitational constant [cm^3/g/s^2]
e      = 4.8032068e-10 # electron charge [esu]
m_e    = 9.1093897e-28 # electron mass [g]
m_p    = 1.6726231e-24 # proton mass [g]
m_n    = 1.6749286-24  # neutron mass [g]
m_H    = 1.660539e-24  # hydrogen mass [g]
amu    = 1.6605402e-24 # atomic mass unit [g]
N_A    = 6.0221367e23  # avagadro's number
k_B    = 1.3806490e-16 # boltzmann constant [erg/K]
eV     = 1.6021766e-12 # electron volt [erg]
keV    = 1e3*eV        # kiloelectron volt [erg]
MeV    = 1e6*eV        # megaelectron volt [erg]
a_rad  = 7.5657233e-15 # radiation density constant [erg/cm^3/K^4]
sig_SB = 5.67051e-5    # stefan-boltzmann constant [erg/cm^2/K^4/s]
alpha  = 7.29735308e-3 # fine-structure constant
Ry     = 2.1798741e-11 # rydberg constant [erg]
sig_T  = 6.6524587e-25 # Thompson scattering cross section [cm^2]

# other useful quantities
M_sol = 1.9891e+33 # solar mass [g]
R_sol = 6.96e10    # solar radius [cm]
L_sol = 3.828e+33  # solar luminosity [erg/s]
T_sol = 5.780e3    # solar temperature [L]
X_sol, Y_sol, Z_sol = 0.7381, 0.2485, 0.0134 # solar abundances
