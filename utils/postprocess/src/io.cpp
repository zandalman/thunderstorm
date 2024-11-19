// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdio>
#include <array>
#include <mpi.h>

// headers
#include "const.h"
#include "functions.h"
#include "parser.h"
#include "io.h"
#include "process.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;
template <typename T>
using vector3d = std::vector<vector2d<T>>;

/// @brief A constructor to initialize the Stat structure.
Stat::Stat(size_t size_, std::string name_, std::string description_)
  : size(size_)
  , name(name_)
  , description(description_)
{}

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
 * @param oss       The string stream.
 * @param vec       The vector to write.
 * @param idx_start The index to start writing the vector.
 * @param idx_end   The index to stop writing the vector.
 */
template <typename T>
void writeVector(
  std::ostringstream &oss, 
  const std::vector<T>& vec, 
  size_t idx_start, 
  size_t idx_end
) {
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
 * @param bin_list      The list of bins.
 * @param stat_list     The list of statistics.
 */
void writeInfo(
  const std::string &infofile_name, 
  Config &config, 
  const vector2d<double> &bin_list,
  const std::vector<Stat> &stat_list
) {

  std::ostringstream oss;

  oss << "Thunderstorm post-processing" << std::endl;
  oss << std::endl;

  oss << "Post-processing parameters" << std::endl;
  oss << "File number:              " << config["IO"]["num_file"] << std::endl;
  oss << "Events per chunk:         " << config["IO"]["num_event_per_chunk"] << std::endl;
  oss << "Number of histories:      " << config["IO"]["num_hist"] << std::endl;
  oss << "Number of densities:      " << config["Grid.Rho"]["num"] << std::endl;
  oss << "Number of Mach numbers:   " << config["Grid.Mach"]["num"] << std::endl;
  oss << "Number of energies:       " << config["Grid.Ener"]["num"] << std::endl;
  oss << "Number of escape lengths: " << config["Grid.Escape"]["num"] << std::endl;
  oss << "Number of lines:          " << 4 + 4 * stat_list.size() << std::endl;
  oss << "Minimum energy [eV]:      " << config["Parameters"]["ener_min"] << std::endl;
  oss << "Turbulence scale [cm]:    " << config["Parameters"]["L"] << std::endl;
  oss << "Geometry:                 " << config["Parameters"]["geo"] << std::endl;
  oss << std::endl;

  oss << "Variable lists" << std::endl;
  oss << "Density" << std::endl;  
  writeVector(oss, bin_list[bin_tag::rho]);
  oss << "Alfven Mach number" << std::endl;  
  writeVector(oss, bin_list[bin_tag::mach]);
  oss << "Energy [eV]" << std::endl;  
  writeVector(oss, bin_list[bin_tag::ener]);
  oss << "Escape [cm]" << std::endl;
  writeVector(oss, bin_list[bin_tag::escape]);
  oss << "Secondary energy [eV]" << std::endl;  
  writeVector(oss, bin_list[bin_tag::ener_sec]);
  oss << "Time [s]" << std::endl;  
  writeVector(oss, bin_list[bin_tag::time]);
  oss << std::endl;

  size_t num_line;
  std::array<std::string, 4> stattype_short_list = {"mean", "var", "skew", "kurt"};
  std::array<std::string, 4> stattype_list = {"mean", "variance", "skewness", "kurtosis"};

  oss << "Post-processed data names" << std::endl;
  oss << "1.  rho" << std::endl;
  oss << "2.  mach_A" << std::endl;
  oss << "3.  ener" << std::endl;
  oss << "4.  escape" << std::endl;
  num_line = 5;
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    for ( size_t j = 0; j < 4; j++ ) {
      std::string space = num_line < 10 ? ".  " : ". ";
      oss << num_line + 0 << space << stat_list[i].name << "_" << stattype_short_list[j] << std::endl;
      num_line++;
    }
  }
  oss << std::endl;

  oss << "Post-processed data descriptions" << std::endl;
  oss << "1.  Density" << std::endl;
  oss << "2.  Alfven Mach number" << std::endl;
  oss << "3.  energy [ev]" << std::endl;
  oss << "4.  escape length [cm]" << std::endl;
  num_line = 5;

  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    for ( size_t j = 0; j < 4; j++ ) {
      std::string space = num_line < 10 ? ".  " : ". ";
      oss << num_line + 0 << space << stattype_list[j] << " " << stat_list[i].description << std::endl;
      num_line++;
    }
  }

  std::ofstream infofile(infofile_name);
  writeOss(infofile, oss);
}

/**
 * @brief Write processed data to an output file.
 * 
 * @param outfile_name        The name of the output file.
 * @param bin_list            The list of bins.
 * @param stat_list           The list of statistics.
 * @param avg_stat_list_flat  The flattened list of mean statistics.
 * @param var_stat_list_flat  The flattened list of variance statistics.
 * @param skew_stat_list_flat The flattened list of skewness statistics.
 * @param kurt_stat_list_flat The flattened list of kurtosis statistics.
 */
void writeData(
  const std::string &outfile_name,
  const vector2d<double> &bin_list, 
  const std::vector<Stat> &stat_list,
  const std::vector<double> &mean_stat_list_flat, 
  const std::vector<double> &var_stat_list_flat,
  const std::vector<double> &skew_stat_list_flat,
  const std::vector<double> &kurt_stat_list_flat
) {
  std::ostringstream oss; // data stream
  size_t idx = 0; // index in flat data vectors

  for ( size_t i = 0; i < bin_list[bin_tag::mach].size(); i++ ) {
    for ( size_t j = 0; j < bin_list[bin_tag::ener].size(); j++ ) {
      for ( size_t k = 0; k < bin_list[bin_tag::escape].size(); k++ ) {
        oss << bin_list[bin_tag::mach][i] << std::endl;
        oss << bin_list[bin_tag::ener][j] << std::endl;
        oss << bin_list[bin_tag::escape][k] << std::endl;
        for ( size_t l = 0; l < stat_list.size(); l++ ) {
          const Stat &stat = stat_list[l];
          writeVector(oss, mean_stat_list_flat, idx, idx + stat.size);
          writeVector(oss, var_stat_list_flat, idx, idx + stat.size);
          writeVector(oss, skew_stat_list_flat, idx, idx + stat.size);
          writeVector(oss, kurt_stat_list_flat, idx, idx + stat.size);
          idx += stat.size; // increment the index
        }
      }
    }
  }

  std::ofstream outfile(outfile_name, std::ios::app);
  writeOss(outfile, oss);
}

/**
 * @brief Write the particle history to an output file.
 * 
 * @param histfile_name The name of the output file.
 * @param bin_list      The list of bins.
 * @param data_grid     The grid of data.
 */
void writeHist(
  const std::string &histfile_name,
  const vector2d<double> &bin_list,
  vector3d<Data>& data_grid
) {
  
  std::ostringstream oss; // data stream
  oss << "format" << std::endl;
  oss << "start: idx_mach, idx_ener, idx_escape" << std::endl;
  oss << "mach, ener[eV], escape[cm]" << std::endl;
  oss << "time[s], x[cm], y[cm], z[cm], cos_alpha, ener[eV], flag" << std::endl;
  oss << "end" << std::endl << std::endl;

  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
    
        oss << "start: " << i << ", " << j << ", " << k << std::endl;
        oss << bin_list[bin_tag::mach][i] << ", ";
        oss << bin_list[bin_tag::ener][j] << ", ";
        oss << bin_list[bin_tag::escape][k] << std::endl;
        oss << data_grid[i][j][k].oss.str();
        oss << "end" << std::endl;

      }
    }
  }

  std::ofstream histfile(histfile_name);
  writeOss(histfile, oss);
}

/** 
 * @brief Write a data stream to a file.
 * 
 * @param outfile The outfile.
 * @param oss     The data stream.
 */
void writeOss(
  std::ofstream &outfile, 
  const std::ostringstream &oss
) {
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
