#ifndef PROCESS_H
#define PROCESS_H

// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <chrono>
#include <memory>

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
  double lam_turb;
  double dt;                       // The zone timestep [s].
  bool super;                      // Whether the turbulence is super-Alfvenic
  size_t ndim;                     // The number of dimensions.
  size_t nker;                     // The number of kernels.
  size_t nmom;                     // The number of moments.
  size_t nstat;                    // The number of statistics.

  double ell_A;
  double gam_par;
  double cut_par;
  double gam_perp;
  double cut_perp;

  double ener;                     // The initial energy [eV].
  double ener_prev;                // The previous energy in Lightning [eV].
  double ener_start;               // The initial energy in Lightning [eV].
  
  double time;                     // The initial time [s].
  double time_prev;                // The previous time in Lightning [s].
  double time_start;               // The initial time in Lightning [s].
  
  bool outoftime;                  // Whether the particle is out of time.
  bool thermalized;                // Whether the particle thermalized.
  
  double s_start;
  double rpar;
  
  Vec pos;                         // The position [cm].
  double splus_prev;               // The previous positive distance along the field line [cm].
  double sminus_prev;              // The previous negative distance along the field line [cm].
  double s_scat;                   // The distance along a field line to scattering [cm].
  
  std::ostringstream oss;          // The string stream.
  vector3d<double> part_stat_list; // The statistics for a single particle.
  vector3d<double> M1_stat_list;   // The mean (M1) statistics.
  vector3d<double> M2_stat_list;   // The M2 statistics.
  vector3d<double> M3_stat_list;   // The M3 statistics.
  vector3d<double> M4_stat_list;   // The M4 statistics.

  Data() = default;
  
  // Move constructor and move assignment
  Data(Data&&) = default;
  Data& operator=(Data&&) = default;

  // Delete copy constructor and copy assignment
  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;

  Data(
    double mach_A_, 
    double dx_, 
    double ener_low_,
    double ener_high_,
    double ener_min_,
    double lam_turb_,
    double dt_,
    int ndim_,
    int nmom_,
    const std::vector<Stat> &stat_list
  );
  void reset();
  void calcStat(int n_int, const std::vector<Stat> &stat_list);
};

void postProcPart(
  int count,
  const vector2d<double> &bin_list,
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
  const size_t nker,
  const size_t nmom,
  std::vector<double> &M1_stat_list_flat,
  std::vector<double> &M2_stat_list_flat, 
  std::vector<double> &M3_stat_list_flat, 
  std::vector<double> &M4_stat_list_flat, 
  size_t &size_flat
);

#endif
