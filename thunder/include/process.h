#ifndef PROCESS_H
#define PROCESS_H

// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <chrono>

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
  double dx;                       // The zone size [cm].
  double ener_low;                 // The lower bound of the energy bin [eV].
  double ener_high;                // The upper bound of the energy bin [eV].
  double ener_min;                 // The minimum energy [eV].
  double dt;                       // The zone timestep [s].
  bool super;                      // Whether the turbulence is super-Alfvenic
  int ndim;                        // The number of dimensions.

  double ener;                     // The current energy [eV].
  double ener_prev;                // The previous energy [eV].
  double ener_start;               // The start energy [eV].
  
  double time;                     // The curret time [s].
  double time_prev;                // The previous time [s].
  double time_start;               // The start time [s].
  bool escaped;                    // Whether the particle has escaped.
  
  Vec pos;                         // The position [cm].
  double splus_prev;               // The previous positive distance along the field line [cm].
  double sminus_prev;              // The previous negative distance along the field line [cm].
  
  Vec Bhat;                        // The B-field direction.
  double lam_scat;                 // The mean free path along a field line to scatter [cm].
  double s_scat;                   // The distance along a field line to scattering [cm].
  
  std::ostringstream oss;          // The string stream.
  vector2d<double> part_stat_list; // The statistics for a single particle.
  vector2d<double> mean_stat_list; // The mean (M1) statistics.
  vector2d<double> M2_stat_list;   // The M2 statistics.
  vector2d<double> M3_stat_list;   // The M3 statistics.
  vector2d<double> M4_stat_list;   // The M4 statistics.

  Data() = default;
  Data(
    double mach_A_, 
    double dx_, 
    double ener_low_,
    double ener_high_,
    double ener_min_,
    double dt_,
    int ndim_,
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
  std::chrono::steady_clock::time_point start,
  int walltime,
  const vector2d<double> &bin_list, 
  const std::vector<Stat>& stat_list, 
  const std::string &histdir_name,
  int idx_hist_max,
  int &idx_hist,
  int &count, 
  vector3d<Data>& data_grid,
  bool &no_time
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
