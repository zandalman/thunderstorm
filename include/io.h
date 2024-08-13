#ifndef IO_H
#define IO_H

// includes
#include <string>

// headers
#include "parser.h"

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
  Event(int id_, int nstep_);
};

void clearInfo(const std::string& infofile_name);
void clearOutfile(const std::string& outfile_name);
void writeInfo(const std::string& infofile_name, int size, Config& config, const Vector1d& ab, const EEDLData& eedl, double n_i, double n_e_free, double lam_deb, double B0);
void writeEvent(const std::string& outfile_name, std::vector<Event> event_list);

#endif
