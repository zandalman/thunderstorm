#ifndef PROCESS_H
#define PROCESS_H

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A structure to represent particle data.
struct Data {
  double ener;                        // The particle energy [eV].
  double escape;                      // The escape legnth [cm].
  double mach_A;                      // The Alfven Mach number.
  double L;                           // The turbulence injection scale [cm].
  int geo;                            // The geometry tag.
  double time_start;                  // The start time [s].
  double splus_start;                 // The start positive distance along the field line [cm].
  double sminus_start;                // The start negative distance along the field line [cm].
  double time;                        // The time [s].
  double splus;                       // The positive distance along the field line [cm].
  double sminus;                      // The negative distance along the field line [cm].
  double rpar;                        // The mean parallel displacement [cm].
  double varpar;                      // The standard deviation of parallel displacement [cm].
  double varperp;                     // The standard deviation of perpendicular displacement [cm].
  double fesc_min;                    // The minimum survival fraction.
  std::vector<double> ener_sec_list;  // The secondary energy bins [eV].
  std::vector<double> ener_loss_mech; // The energy loss for each mechanism [eV].
  std::vector<double> num_ion_elem;   // The number of ionizations per element.
  std::vector<double> num_sec_ener;   // The number of secondary particles per energy bin.

  Data() = default;
  Data(double ener_, double escape_, double mach_A_, double L_, int geo_, std::vector<double> ener_sec_list_);
  void reset();
};

void processEvent(const Event* event, vector2d<Data>& data_grid);
void processFile(const std::string &datafile_name, size_t num_event_per_chunk, int &count, vector2d<Data>& data_grid);
void getFlatData(const vector2d<Data>& data_grid, std::vector<double> &ener_loss_mech_flat, std::vector<double> &num_ion_elem_flat, std::vector<double> &num_sec_ener_flat);

#endif
