#ifndef SIM_H
#define SIM_H

// includes
#include <string>
#include <vector>

// headers
#include "io.h"
#include "parser.h"
#include "part.h"

/// @brief A structure to represent the simulation.
struct Sim {
  Part part;                     // The particle object.
  const EEDLData& eedl;          // Data from the EEDL database.
  const Vector1d& ab;            // A vector of elemental abundances.
  std::string outfile;           // The outfile name.
  double rho;                    // The density [g/cc].
  double Bmag_co;                // The amplitude of the coherent magnetic field [G].
  double Bmag_turb;              // The amplitude of the turbulent magnetic field [G].
  bool do_Bfield;                // Whether a magnetic field is present.
  int nstep;                     // The step number.
  double time;                   // The time.
  std::vector<Event> event_list; // A vector of event objects.

  Sim(Part part_, const EEDLData& eedl_, const Vector1d& ab_, std::string outfile_, double rho_, double Bmag_co_, double Bmag_turb_);
  double calcSigTot();
  void move(double sig_tot, Event &event);
  int choseElem();
  int choseInter(int Zelem);
  int choseIon(int Zelem);
  void interact(Event &event);
  void step();
  void multistep(int num);
};

#endif
