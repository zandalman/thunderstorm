#ifndef IO_H
#define IO_H

// includes
#include <vector>
#include <string>

// headers
#include "parser.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A structure to represent a statistics.
struct Stat {
  size_t size;
  std::string name;
  std::string description;

  Stat() = default;
  Stat(size_t size_, std::string name_, std::string description_);
};

/// @brief A structure to represent an event.
struct Event {
  int id;                  // The particle ID.
  int nstep;               // The step number.
  int Zelem;               // The proton number of the element.
  int interaction;         // The interaction flag.
  int ion;                 // The ion index.
  double time;             // The event time [s].
  double splus;            // The positive distance along the field line [cm].
  double sminus;           // The negative distance along the field line [cm].
  double ener;             // The particle kinetic energy [eV].
  double cos_alpha;        // The particle pitch angle cosine.
  double cos_th;           // The cosine of the scattering angle.
  double ener_loss;        // The energy lost [eV].
  double ener_sec;         // The energy of the secondary [eV].
  double ener_loss_sync;   // The energy lost due to synchrotron [eV].
  double ener_loss_cher;   // The energy lost due to Cherenkov radiation [eV].
  double ener_loss_moller; // The energy lost due to small-angle Moller scattering [eV].

  Event() = default;
  Event(int int_data[5], double double_data[11]);
};

template <typename T>
void writeVector(
  std::ostringstream &oss, 
  const std::vector<T>& vec, 
  size_t idx_start = 0, 
  size_t idx_end = 0
);
void clearFile(const std::string& file_name);
void writeInfo(
  const std::string &infofile_name, 
  Config &config, 
  const vector2d<double> &bin_list,
  const std::vector<Stat> &stat_list
);
void writeData(
  const std::string &outfile_name,
  const vector2d<double> &bin_list, 
  const std::vector<Stat> &stat_list,
  const std::vector<double> &mean_stat_list_flat, 
  const std::vector<double> &var_stat_list_flat,
  const std::vector<double> &skew_stat_list_flat,
  const std::vector<double> &kurt_stat_list_flat
);

#endif
