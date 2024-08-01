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

/// @brief A structure to represent an event.
struct Event {
  int id;                  // The particle ID.
  int nstep;               // The step number.
  int Zelem;               // The proton number of the element.
  int interaction;         // The interaction flag.
  int ion;                 // The ion index.
  double time;             // The event time [s].
  double x, y, z;          // The event coordinates [cm].
  double ener;             // The particle kinetic energy [eV].
  double cos_th;           // The cosine of the scattering angle.
  double cos_alpha;        // The particle pitch angle cosine.
  double ener_loss;        // The energy lost [eV].
  double ener_sec;         // The energy of the secondary [eV].
  double ener_loss_sync;   // The energy lost due to synchrotron [eV].
  double ener_loss_cher;   // The energy lost due to Cherenkov radiation [eV].
  double ener_loss_moller; // The energy lost due to small-angle Moller scattering [eV].

  Event() = default;
  Event(int int_data[5], double double_data[12]);
};

/// @brief A structure to represent post-processed particle data.
struct PartData {
  int id;                             // The particle ID.
  int num_ev = 0;                     // The number of events.
  double ener = 0.;                   // The particle energy [eV].
  double t_start = 0.;                // The start time [s].
  double x_start = 0.;                // The start x-coordinate [cm].
  double y_start = 0.;                // The start y-coordinate [cm].
  double z_start = 0.;                // The start z-coordinate [cm].
  double t_end = 0.;                  // The death time [s].
  double x_end = 0.;                  // The death x-coordinate [cm].
  double y_end = 0.;                  // The death y-coordinate [cm].
  double z_end = 0.;                  // The death z-coordinate [cm].
  std::vector<double> ener_loss_mech; // The energy loss for each mechanism [eV].
  std::vector<int> num_ion_elem;      // The number of ionizations per element.
  std::vector<int> num_sec_ener;      // The number of secondary particles per energy bin.
  std::vector<int> num_ev_time;       // The number of events per time bin.
  std::vector<double> dis_par_time;   // The average parallel distance in each time bin.
  std::vector<double> dis_perp_time;  // The average perpdicular distance in each time bin.
  std::vector<double> ener_time;      // The average energy in each time bin.
  std::vector<double> ener_loss_time; // The energy loss per time bin.
  vector2d<double> ener_loss_dis2d;   // The energy loss per parallel and perpendicular distance bin.

  PartData() = default;
  PartData(size_t num_ener_sec, size_t num_time, size_t num_dis_par, size_t num_dis_perp, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
  void reset(int id_, double ener_, size_t num_ener_sec, size_t num_time, size_t num_dis, size_t num_dis_perp, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
};

template <typename T>
void writeVector(std::ostringstream &oss, const std::vector<T>& vec);
template <typename T>
void writeVector2d(std::ostringstream &oss, const vector2d<T>& vec);
void clearFile(const std::string& file_name);
void writeInfo(std::string &infofile_name, Config& config, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list, size_t num_stat, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
void processEvent(const Event* event, std::vector<PartData>& part_data_list, const std::vector<double>& ener_list, const std::vector<double> &ener_sec_list, const std::vector<double>& time_list, const std::vector<double>& dis_par_list, const std::vector<double>& dis_perp_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
void postProcPartData(std::vector<PartData>& part_data_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
void writePartData(std::string &outfile_name, std::vector<PartData>& part_data_list, const std::vector<bool> &do_stat_list);
int processFile(std::string &datafile_name, std::string &outfile_name, size_t num_event_per_chunk, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat);
void concatenateFiles(std::string data_path, int num_file);

#endif
