#ifndef PROCESS_H
#define PROCESS_H

// includes
#include <vector>

// headers
#include "io.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A structure to represent data.
struct Data {
  double ener;                             // The particle energy [eV].
  double escape;                           // The escape legnth [cm].
  double mach_A;                           // The Alfven Mach number.
  double L;                                // The turbulence injection scale [cm].
  int geo;                                 // The geometry tag.
  double time_start;                       // The start time [s].
  double splus_start;                      // The start positive distance along the field line [cm].
  double sminus_start;                     // The start negative distance along the field line [cm].
  double time;                             // The time [s].
  double splus;                            // The positive distance along the field line [cm].
  double sminus;                           // The negative distance along the field line [cm].
  double rpar;                             // The mean parallel displacement [cm].
  double varpar;                           // The standard deviation of parallel displacement [cm].
  double varperp;                          // The standard deviation of perpendicular displacement [cm].
  double fesc_min;                         // The minimum survival fraction.
  vector2d<double> part_stat_list;
  vector2d<double> mean_stat_list;
  vector2d<double> M2_stat_list;
  vector2d<double> M3_stat_list;
  vector2d<double> M4_stat_list;

  Data() = default;
  Data(double ener_, double escape_, double mach_A_, double L_, int geo_, const std::vector<Stat> &stat_list);
  void reset();
  void calcStat(int n_int, const std::vector<Stat> &stat_list);
};

void processEvent(const Event* event, const vector2d<double> &bin_list, const std::vector<Stat>& stat_list, vector2d<Data>& data_grid);
void processFile(const std::string &datafile_name, size_t num_event_per_chunk, const vector2d<double> &bin_list, const std::vector<Stat>& stat_list, int &count, vector2d<Data>& data_grid);
void getFlatData(const vector2d<Data>& data_grid, const std::vector<Stat> &stat_list, std::vector<double> &mean_stat_list_flat, std::vector<double> &M2_stat_list_flat, std::vector<double> &M3_stat_list_flat, std::vector<double> &M4_stat_list_flat, size_t &size_flat);

#endif
