#ifndef IO_H
#define IO_H

// includes
#include <string>

/// @brief A structure to represent an event.
struct Event {
  int id;                // The particle ID.
  int nstep;             // The step number.
  int Zelem;             // The proton number of the element.
  int interaction;       // The interaction flag.
  int ion;               // The ion index.
  double time;           // The event time.
  double x, y, z;        // The event coordinates.
  double ener;           // The particle energy.
  double cos_th;         // The cosine of the scattering angle.
  double ener_loss;      // The energy lost.
  double ener_sec;       // The energy of the secondary.
  double ener_loss_sync; // The energy lost due to synchrotron.

  Event() = default;
  Event(int id_, int nstep_);
};

void clearOutfile(const std::string& outfile);
void writeEvent(std::string outfile, std::vector<Event> event_list);

#endif