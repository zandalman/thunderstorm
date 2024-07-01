#ifndef IO_H
#define IO_H

// includes
#include <vector>
#include <string>

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
  double ener_loss;        // The energy lost [eV].
  double ener_sec;         // The energy of the secondary [eV].
  double ener_loss_sync;   // The energy lost due to synchrotron [eV].
  double ener_loss_cher;   // The energy lost due to Cherenkov radiation [eV].
  double ener_loss_moller; // The energy lost due to small-angle Moller scattering [eV].

  Event() = default;
  Event(int int_data[5], double double_data[11]);
};

struct PartData {
  int id;
  int num_ev = 0;
  double ener = 0.;
  double x_start = 0.;
  double y_start = 0.;
  double z_start = 0.;
  double ener_loss_mech[6];
  int num_ion[118];
  std::vector<int> num_ev_time;
  std::vector<double> ener_time;
  std::vector<double> ener_loss_time;
  std::vector<double> ener_loss_dis;

  PartData() = default;
  PartData(size_t num_time, size_t num_dis);
  void reset(int id_, double ener_, size_t num_time, size_t num_dis);
};

void processEvent(const Event* event, std::vector<PartData>& part_data_list, const std::vector<double>& ener_list, const std::vector<double>& time_list, const std::vector<double>& dis_list);
void writePartData(std::vector<PartData>& part_data_list, const std::vector<double> &ener_list, const std::vector<double> &time_list, const std::vector<double> &dis_list);
int processFile(std::string filename, size_t num_event_per_chunk, const std::vector<double> &ener_list, const std::vector<double> &time_list, const std::vector<double> &dis_list);

#endif
