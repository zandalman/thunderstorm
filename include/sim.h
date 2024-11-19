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
  double temp;                   // The temperature [K].
  double ion_state_avg;          // The average ionization state.
  double cos_th_cut;             // The cutoff scattering angle cosine for discrete Moller scattering.
  int nstep;                     // The step number.
  double time;                   // The simulation time [s].
  double n_i;                    // The ion number density [1/cc].
  double n_e_free;               // The free electron number density [1/cc].
  double lam_deb;                // The Debye length [cm].
  double B0;                     // The coherent magnetic field amplitude [G].
  bool do_moller;                // Do moller scattering and energy losses.
  bool do_cerenkov;              // Do Cerenkov energy losses.
  bool do_sync;                  // Do synchrotron energy losses.
  std::vector<Event> event_list; // A vector of event objects.

  Sim(
    Part part_, 
    const EEDLData& eedl_, 
    const Vector1d& ab_, 
    std::string outfile_, 
    double rho_, 
    double temp_, 
    double ion_state_avg_, 
    double B0_, 
    double cos_th_cut_,
    bool do_moller_,
    bool do_cerenkov_,
    bool do_sync_
  );
  void reset(Part part);
  void kill();
  double calcSigTot();
  void move(double sig_tot, Event &event);
  int choseElem();
  int choseInter(int Zelem);
  int choseIon(int Zelem);
  void interact(Event &event);
  void step();
};

#endif
