// includes
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>

// headers
#include "io.h"
#include "parser.h"

/// @brief A constructor to initialize the Event structure.
Event::Event(int id_, int nstep_) 
  : id(id_)                // The particle ID.
  , nstep(nstep_)          // The step number.
  , Zelem(-1)              // The proton number of the element.
  , interaction(-1)        // The interaction flag.
  , ion(-1)                // The ion index.
  , time(0.0)              // The event time [s].
  , x(0.0), y(0.0), z(0.0) // The event coordinates [cm].
  , ener(0.0)              // The particle kinetic energy [eV].
  , cos_alpha(1.0)         // The particle pitch angle cosine.
  , cos_th(1.0)            // The cosine of the scattering angle.
  , ener_loss(0.0)         // The energy lost [eV].
  , ener_sec(0.0)          // The energy of the secondary [eV].
  , ener_loss_sync(0.0)    // The energy lost due to synchrotron [eV].
  , ener_loss_cher(0.0)    // The energy lost due to Cherenkov radiation [eV].
  , ener_loss_moller(0.0)  // The energy lost due to small-angle Moller scattering [eV].
 {}

/**
 * @brief Clear the info file.
 * 
 * @param outfile The info file name.
*/
void clearInfo(const std::string& infofile_name) {
  std::ofstream infofile(infofile_name);
  infofile.close();
}

/**
 * @brief Clear the outfile.
 * 
 * @param outfile The outfile name.
*/
void clearOutfile(const std::string& outfile_name) {
    std::ofstream outfile(outfile_name, std::ios::binary | std::ios::trunc);
    outfile.close();
}

/**
 * @brief Write information about the simulation to a text file.
 * 
 * @param outfile The info file name.
*/
void writeInfo(const std::string& infofile_name, int size, Config& config, const Vector1d& ab, const EEDLData& eedl, double n_i, double n_e_free, double lam_deb, double B0) {

  std::ofstream infofile(infofile_name);
  std::ifstream licensefile("../LICENSE");

  if ( !infofile ) {
    std::cerr << "Failed to open " << infofile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if ( !licensefile ) {
    std::cerr << "Failed to open " << "../LICENSE" << " for reading." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  infofile << "Thunderstorm" << std::endl;
  std::string line;
  while ( std::getline(licensefile, line) ) {
    infofile << line << std::endl;
  }
  infofile << std::endl;

  infofile << "Simulation parameters" << std::endl;
  infofile << "Number of MPI processes:         " << size << std::endl;
  infofile << "Simulation duration [s]:         " << config["Simulation"]["tmax"] << std::endl;
  infofile << "Density [g/cc]:                  " << config["Background"]["rho"] << std::endl;
  infofile << "Temperature [K]:                 " << config["Background"]["temp"] << std::endl;
  infofile << "Abundance time [day]:            " << config["Background"]["ab_time"] << std::endl;
  infofile << "Average ion state:               " << config["Background"]["ion_state_avg"] << std::endl;
  infofile << "Ion number density [1/cc]:       " << n_i << std::endl;
  infofile << "Free elec number density [1/cc]: " << n_e_free << std::endl;
  infofile << "Debye length [cm]:               " << lam_deb << std::endl;
  infofile << "Particle energy [eV]:            " << config["Particle"]["ener"] << std::endl;
  infofile << "Particle lifetime [s]:           " << config["Particle"]["tpart"] << std::endl;
  infofile << "Turbulence injection scale [cm]: " << config["Bfield"]["L"] << std::endl;
  infofile << "Plasma beta:                     " << config["Bfield"]["beta"] << std::endl;
  infofile << "Coherent B-field amplitude [G]:  " << B0 << std::endl;
  infofile << "Alfven Mach number:              " << config["Bfield"]["mach_A"] << std::endl;
  infofile << "Turbulent cross section fraction " << config["Bfield"]["sig_turb_frac"] << std::endl;
  infofile << "B-field curvature spectrum expon " << config["Bfield"]["alpha"] << std::endl;
  infofile << "Discrete Moller cos angle cutoff " << config["Simulation"]["cos_th_cut"] << std::endl;
  infofile << std::endl;

  infofile << "Abundances" << std::endl;
  infofile << "Z,ab" << std::endl;
  for ( size_t i = 0; i < ab.size(); i++ ) {
    infofile << i << "," << ab[i] << std::endl;
  }
  infofile << std::endl;

  infofile << "Ions" << std::endl;
  infofile << "Z:ion_list" << std::endl;
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    infofile << i+1 << ":";
    std::vector<std::string> ion_list = eedl[i].ion_list;
    for ( size_t j = 0; j < ion_list.size(); j++ ) {
      infofile << ion_list[j] << ",";
    }
    infofile << std::endl;
  }
  infofile << std::endl;

  if ( !infofile.good() ) {
    std::cerr << "Error writing to " << infofile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  licensefile.close();
  infofile.close();
}

/**
 * @brief Write a list of events to the outfile.
 * 
 * @param outfile    The outfile name.
 * @param event_list A vector of events.
*/
void writeEvent(const std::string& outfile_name, std::vector<Event> event_list) {
  std::ofstream outfile(outfile_name, std::ios::out | std::ios::binary | std::ios::app);

  if ( !outfile ) {
    std::cerr << "Failed to open " << outfile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  outfile.write(reinterpret_cast<const char*>(event_list.data()), event_list.size() * sizeof(Event));
  
  if ( !outfile.good() ) {
    std::cerr << "Error writing to " << outfile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  outfile.close();
}
