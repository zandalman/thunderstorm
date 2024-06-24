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
  , time(0.0)              // The event time.
  , x(0.0), y(0.0), z(0.0) // The event coordinates.
  , ener(0.0)              // The particle energy.
  , cos_th(1.0)            // The cosine of the scattering angle.
  , ener_loss(0.0)         // The energy lost.
  , ener_sec(0.0)          // The energy of the secondary.
  , ener_loss_sync(0.0)    // The energy lost due to synchrotron.
 {}

/**
 * @brief Clear the info file.
 * 
 * @param outfile The info file name.
*/
void clearInfo(const std::string& infofile) {
  std::ofstream file(infofile);
  file.close();
}

/**
 * @brief Clear the outfile.
 * 
 * @param outfile The outfile name.
*/
void clearOutfile(const std::string& outfile) {
    std::ofstream file(outfile, std::ios::binary | std::ios::trunc);
    file.close();
}

/**
 * @brief Write information about the simulation to a text file.
 * 
 * @param outfile The info file name.
*/
void writeInfo(const std::string& infofile, Config& config, const Vector1d& ab, const EEDLData& eedl) {

  std::ofstream file;
  std::ifstream licenseFile("../LICENSE");
  file.open(infofile);

  if ( !file ) {
    std::cerr << "Failed to open info file for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if ( !licenseFile ) {
    std::cerr << "Failed to open license file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  file << "Thunderstorm" << std::endl;
  std::string line;
  while ( std::getline(licenseFile, line) ) file << line << std::endl;
  file << std::endl;

  file << "Simulation parameters" << std::endl;
  file << "Number of particles:         " << config["Simulation"]["npac"] << std::endl;
  file << "Max simulation duration:     " << config["Simulation"]["tmax"] << " [s]" << std::endl;
  file << "Density:                     " << config["Simulation"]["rho"] << " [g/cc]" << std::endl;
  file << "Abundance time:              " << config["Simulation"]["ab_time"] << " [day]" << std::endl;
  file << "Particle energy:             " << config["Particle"]["ener"] << " [eV]" << std::endl;
  file << "Particle lifetime:           " << config["Particle"]["tpart"] << " [s]" << std::endl;
  file << "Coherent B_field amplitude:  " << config["Bfield"]["Bmag_co"] << " [G]" << std::endl;
  file << "Turbulent B-field amplitude: " << config["Bfield"]["Bmag_turb"] << " [G]" << std::endl;
  file << "B-field spectrum index:      " << config["Bfield"]["q"] << std::endl;
  file << "B-field spectrum max scale:  " << config["Bfield"]["Lmax"] << " [cm]" << std::endl;
  file << std::endl;

  file << "Abundances" << std::endl;
  file << "Z,ab" << std::endl;
  for ( size_t i = 0; i < ab.size(); i++ ) {
    file << i << "," << ab[i] << std::endl;
  }
  file << std::endl;

  file << "Ions" << std::endl;
  file << "Z:ion_list" << std::endl;
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    file << i+1 << ":";
    std::vector<std::string> ion_list = eedl[i].ion_list;
    for ( size_t j = 0; j < ion_list.size(); j++ ) {
      file << ion_list[j] << ",";
    }
    file << std::endl;
  }
  file << std::endl;

  licenseFile.close();
  file.close();

  if ( !file.good() ) {
    std::cerr << "Error writing to info file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if ( !licenseFile.good() ) {
    std::cerr << "Error closing license file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

/**
 * @brief Write a list of events to the outfile.
 * 
 * @param outfile    The outfile name.
 * @param event_list A vector of events.
*/
void writeEvent(std::string outfile, std::vector<Event> event_list) {
  std::ofstream file(outfile, std::ios::out | std::ios::binary | std::ios::app);

  if ( !file ) {
    std::cerr << "Failed to open file for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  file.write(reinterpret_cast<const char*>(event_list.data()), event_list.size() * sizeof(Event));
  file.close();

  if ( !file.good() ) {
    std::cerr << "Error writing to file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}
