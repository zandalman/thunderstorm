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
  Vector1d ion_state;            // The ionization state vector.
  double q_avg;                  // The average ionization state.
  double qsq_avg;                // The average square ionization state.
  double cos_th_cut;             // The cutoff scattering angle cosine for discrete Moller scattering.
  int nstep;                     // The step number.
  double time;                   // The simulation time [s].
  double n_i;                    // The ion number density [1/cc].
  double n_e_free;               // The free electron number density [1/cc].
  double lam_deb;                // The Debye length [cm].
  double B0;                     // The coherent magnetic field amplitude [G].
  bool neutral;                  // Whether the ejecta is neutral.
  bool do_cerenkov;              // Do Cerenkov energy losses.
  bool do_sync;                  // Do synchrotron energy losses.
  double mmw;                    // Mean molecular weight [g/mol].
  std::vector<Event> event_list; // A vector of event objects.

  Sim(
    Part part_, 
    const EEDLData& eedl_, 
    const Vector1d& ab_, 
    std::string outfile_, 
    double rho_, 
    double temp_, 
    Vector1d ion_state_, 
    double B0_, 
    double cos_th_cut_,
    bool neutral_,
    bool do_cerenkov_,
    bool do_sync_,
    double mmw_
  );
  void reset(Part part);
  void kill();
  double calcSigTot();
  void move(double sig_tot, Event &event);
  void choseElem(bool &moller, int &Zelem, int &stage);
  int choseInter(int Zelem, int stage);
  int choseIon(int Zelem, int stage);
  void interact(Event &event);
  void step();
};

#endif
