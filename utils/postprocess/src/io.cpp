// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdio>
#include <mpi.h>

// headers
#include "const.h"
#include "functions.h"
#include "parser.h"
#include "io.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A constructor to initialize the Event structure.
Event::Event(int int_data[5], double double_data[11]) 
  : id(int_data[0])                   // The particle ID.
  , nstep(int_data[1])                // The step number.
  , Zelem(int_data[2])                // The proton number of the element.
  , interaction(int_data[3])          // The interaction flag.
  , ion(int_data[4])                  // The ion index.
  , time(double_data[0])              // The event time [s].
  , splus(double_data[1])             // The positive distance along the field line [cm].
  , sminus(double_data[2])            // The negative distance along the field line [cm].
  , ener(double_data[3])              // The particle kinetic energy [eV].
  , cos_alpha(double_data[4])         // The particle pitch angle cosine.
  , cos_th(double_data[5])            // The cosine of the scattering angle.
  , ener_loss(double_data[6])         // The energy lost [eV].
  , ener_sec(double_data[7])          // The energy of the secondary [eV].
  , ener_loss_sync(double_data[8])    // The energy lost due to synchrotron [eV].
  , ener_loss_cher(double_data[9])    // The energy lost due to Cherenkov radiation [eV].
  , ener_loss_moller(double_data[10]) // The energy lost due to small-angle Moller scattering [eV].
 {}


/**
 * @brief Clear a file.
 * 
 * @param file_name The name of the file to clear.
*/
void clearFile(const std::string& file_name) {
  std::ofstream file(file_name);
  file.close();
}

/**
 * @brief Write a vector to a string stream.
 * 
 * @param oss The string stream.
 * @param vec The vector to write.
 */
template <typename T>
void writeVector(std::ostringstream &oss, const std::vector<T>& vec, size_t idx_start, size_t idx_end) {
  idx_end = idx_end == 0 ? vec.size() : idx_end;
  for (size_t i = idx_start; i < idx_end; i++) {
    oss << vec[i];
    if (i != idx_end - 1) oss << ",";
  }
  oss << std::endl;
}

/**
 * @brief Write information about the post-processing to a text file.
 * 
 * @param infofile_name The name of the information file.
 * @param config        The struct containing the configuration file information.
 * @param ener_list     The list of energy bins [eV].
 * @param ener_sec_list The list of secondary energy bins [eV].
 */
void writeInfo(const std::string &infofile_name, Config &config, const std::vector<double> &ener_list, const std::vector<double> &escape_list, const std::vector<double> &ener_sec_list) {

  std::ostringstream oss;

  oss << "Thunderstorm post-processing" << std::endl;
  oss << std::endl;

  oss << "Post-processing parameters" << std::endl;
  oss << "File number:              " << config["IO"]["num_file"] << std::endl;
  oss << "Events per chunk:         " << config["IO"]["num_event_per_chunk"] << std::endl;
  oss << "Number of energies:       " << config["Bin.Energy"]["num"] << std::endl;
  oss << "Number of escape lengths: " << config["Bin.Escape"]["num"] << std::endl;
  oss << "Number of lines:          " << 5 << std::endl;
  oss << "Alfven Mach number:       " << config["Bfield"]["mach_A"] << std::endl;
  oss << "Turbulence scale [cm]:    " << config["Bfield"]["L"] << std::endl;
  oss << "Geometry:                 " << config["Geometry"]["geo"] << std::endl;
  oss << std::endl;

  oss << "Variable lists" << std::endl;
  oss << "Energy [eV]" << std::endl;  
  writeVector(oss, ener_list);
  oss << "Escape [cm]" << std::endl;
  writeVector(oss, escape_list);
  oss << "Secondary energy [eV]" << std::endl;  
  writeVector(oss, ener_sec_list);
  oss << std::endl;

  size_t num_line = 1;
  oss << "Post-processed data format" << std::endl;
  oss << "1.  energy [ev]" << std::endl;
  oss << "2.  escape length [cm]" << std::endl;
  num_line = 3;

  // list of descriptions of each statistic
  std::vector<std::string> stat_def_list = {
    "energy loss [eV] for each mechanism",
    "number of ionizations per element",
    "number of secondary electrons per secondary energy bin"
  };

  for ( size_t i = 0; i < stat_def_list.size(); i++ ) {
    oss << num_line << ".  " << stat_def_list[i] << std::endl;
    num_line++;
  }

  std::ofstream infofile(infofile_name);
  if ( !infofile ) {
    std::cerr << "Failed to open file " << infofile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  infofile << oss.str();
  infofile.close();
  if ( !infofile.good() ) {
    std::cerr << "Error writing to file " << infofile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

/**
 * @brief Write processed data to an output file.
 * 
 * @param outfile_name        The name of the output file.
 * @param ener_list           The list of energies.
 * @param escape_list         The list of escape lengths.
 * @param num_ener_sec        The number of secondary energies.
 * @param ener_loss_mech_flat The flattened list of energy loss for each mechanism [eV]
 * @param num_ion_elem_flat   The flattened list of ionizations per element.
 * @param num_sec_ener_flat   The flattened list of secondary particles per energy bin.
 */
void writeData(
  const std::string &outfile_name,
  const std::vector<double> &ener_list, const std::vector<double> &escape_list,
  size_t num_ener_sec,
  const std::vector<double> &ener_loss_mech_flat, const std::vector<double> &num_ion_elem_flat, const std::vector<double> &num_sec_ener_flat
) {

  // start the stream
  std::ostringstream oss;

  for ( size_t i = 0; i < ener_list.size(); i++ ) {
    for ( size_t j = 0; j < escape_list.size(); j++ ) {
      oss << ener_list[i] << std::endl;
      oss << escape_list[j] << std::endl;
      writeVector(oss, ener_loss_mech_flat, i * j * num_mech, (i * j + 1) * num_mech);
      writeVector(oss, num_ion_elem_flat, i * j * num_elem, (i * j + 1) * num_elem);
      writeVector(oss, num_sec_ener_flat, i * j * num_ener_sec, (i * j + 1) * num_ener_sec);
    }
  }

  std::ofstream outfile(outfile_name, std::ios::app);
  if ( !outfile ) {
    std::cerr << "Failed to open file for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  outfile << oss.str();
  outfile.close();
  if ( !outfile.good() ) {
    std::cerr << "Error writing to file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}
