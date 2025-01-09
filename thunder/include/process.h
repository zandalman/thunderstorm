#ifndef PROCESS_H
#define PROCESS_H

// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

// headers
#include "io.h"
#include "vec.h"
#include "functions.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;
template <typename T>
using vector3d = std::vector<vector2d<T>>;

/// @brief A structure to represent data.
struct Data {
  double mach_A;                   // The Alfven Mach number.
  double ener;                     // The particle energy [eV].
  double ener_min;                 // The minimum energy [eV].
  double inner;                    // The position of the thermalization barrier [cm].
  double outer;                    // The position of the escape barrier [cm].
  double turb;                     // The turbulence injection scale [cm].
  bool escaped;                    // Whether the particle has escaped.
  bool therm;                      // Whether the particle has thermalized.
  double ener_start;               // The start energy [eV].
  double time_start;               // The start time [s].
  double ener_prev;                // The energy [eV].
  double time_prev;                // The time [s].
  double splus_prev;               // The positive distance along the field line [cm].
  double sminus_prev;              // The negative distance along the field line [cm].
  double lam_scat;                 // The mean free path along a field line to scatter [cm].
  double s_scat;                   // The distance along a field line to scattering [cm].
  Vec pos;                         // The position [cm].
  Vec Bhat;                        // The B-field direction.
  std::ostringstream oss;          // The string stream.
  vector2d<double> part_stat_list; // The statistics for a single particle.
  vector2d<double> mean_stat_list; // The mean statistics
  vector2d<double> M2_stat_list;   // The M2 statistics
  vector2d<double> M3_stat_list;   // The M3 statistics
  vector2d<double> M4_stat_list;   // The M4 statistics

  Data() = default;
  Data(
    double mach_A_, 
    double scale_, 
    double ener_,
    MiscParam misc_param_,
    const std::vector<Stat> &stat_list
  );
  void reset();
  void calcStat(int n_int, const std::vector<Stat> &stat_list);
};

void postProcPart(
  int count,
  const std::vector<Stat> &stat_list, 
  vector3d<Data> &data_grid
);
void processEvent(
  const Event* event, 
  bool do_hist,
  const vector2d<double> &bin_list, 
  vector3d<Data>& data_grid
);
void processFile(
  const std::string &datafile_name, 
  size_t num_event_per_chunk, 
  const vector2d<double> &bin_list, 
  const std::vector<Stat>& stat_list, 
  const std::string &histdir_name,
  int idx_hist_max,
  int &idx_hist,
  int &count, 
  vector3d<Data>& data_grid
);
void getFlatData(
  const vector3d<Data>& data_grid, 
  const std::vector<Stat> &stat_list, 
  std::vector<double> &mean_stat_list_flat,
  std::vector<double> &M2_stat_list_flat, 
  std::vector<double> &M3_stat_list_flat, 
  std::vector<double> &M4_stat_list_flat, 
  size_t &size_flat
);

#endif
