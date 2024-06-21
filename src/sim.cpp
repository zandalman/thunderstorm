// includes
#include <cmath>
#include <algorithm>
#include <iostream>

// headers
#include "const.h"
#include "random.h"
#include "part.h"
#include "vec.h"
#include "sim.h"
#include "parser.h"
#include "functions.h"
#include "io.h"

/// @brief A constructor to initial to the Sim structure.
Sim::Sim(Part part_, const EEDLData& eedl_, const Vector1d& ab_, std::string outfile_, double rho_, double Bmag_co_, double Bmag_turb_, double q_, double Lmax_)
  : part(part_)                               // The particle object.
  , eedl(eedl_)                               // Data from the EEDL database.
  , ab(ab_)                                   // A vector of elemental abundances.
  , outfile(outfile_)                         // The outfile name.
  , rho(rho_)                                 // The density [g/cc].
  , Bmag_co(Bmag_co_)                         // The amplitude of the coherent magnetic field [G].
  , Bmag_turb(Bmag_turb_)                     // The amplitude of the turbulent magneitc field [G].
  , q(q_)                                     // The power law exponent of the magnetic turbulence spectrum.
  , Lmax(Lmax_)                               // The largest scale of magnetic turbulence [cm].
  , do_Bfield(Bmag_turb > 0. || Bmag_co > 0.) // Whether a magnetic field is present.
  , nstep (0)                                 // The step number.
  , time(0.0)                                 // The simulation time [s].
  {
    part.newBvec(Bmag_co, Bmag_turb);
  }

/**
 * @brief Compute the total cross section.
 * 
 * @return The total cross section [cm^2].
*/
double Sim::calcSigTot() {
  double sig_tot = 0.0;
  double sig_Bturb = do_Bfield ? calcSigBturb(part, part.Bvec.mag(), Bmag_turb, rho, q, Lmax) : 0.;
  sig_tot += sig_Bturb;
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    double sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second, true, false, 0., 0.);
    sig_tot += sig * ab[i+1];
  }
  return sig_tot;
}

/**
 * @brief Transport the particle.
 * 
 * @param sig_tot The total cross section [cm^2].
 * @param event   The event structure.
*/
void Sim::move(double sig_tot, Event &event) {
  double dis = -log(1 - xi()) / (rho * constants::N_A * sig_tot);
  double dt = dis / (constants::c * part.beta());
  time += dt;
  if ( do_Bfield ) {
    double cos_alpha = dot(part.vel.unit(), part.Bvec.unit());
    part.pos = part.pos + dis * cos_alpha * part.Bvec.unit();
    part.vel = rotate(part.vel, part.Bvec, cos(2 * M_PI * xi()));
    event.ener_loss_sync = calcPowerSync(part.q_i, part.gam(), part.beta(), cos_alpha, part.Bvec.mag(), part.m_i) * dt;
    part.loseEner(event.ener_loss_sync);
  } else {
    part.pos = part.pos + dis*part.vel.unit();
  }
  event.time = time;
  event.x = part.pos.x;
  event.y = part.pos.y;
  event.z = part.pos.z;
}

/**
 * @brief Select an element.
 * 
 * @return The proton number of the selected element, or 0 for a magnetic field interaction.
*/
int Sim::choseElem() {
  double sig_tot = 0.0;
  Vector1d sig_cum;
  double sig_Bturb = do_Bfield ? calcSigBturb(part, part.Bvec.mag(), Bmag_turb, rho, q, Lmax) : 0.;  
  sig_tot += sig_Bturb;
  sig_cum.push_back(sig_tot);
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    double sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second, true, false, 0., 0.);
    sig_tot += sig * ab[i+1];
    sig_cum.push_back(sig_tot);
  }
  return findIdx(sig_tot*xi(), sig_cum);
}

/**
 * @brief Select an interaction.
 * 
 * @param Zelem The proton number of the element.
 * @return The interaction flag.
*/
int Sim::choseInter(int Zelem) {
  SpecData spec_data = eedl[Zelem-1];
  double sig_tot = 0.0;
  Vector1d sig_cum;
  double sig_scat = interp(part.ener, spec_data.sig_scat_data.first, spec_data.sig_scat_data.second, true, false, 0., 0.);
  double sig_brem = interp(part.ener, spec_data.sig_brem_data.first, spec_data.sig_brem_data.second, true, false, 0., 0.);
  double sig_exc = interp(part.ener, spec_data.sig_exc_data.first, spec_data.sig_exc_data.second, true, false, 0., 0.);
  double sig_ion = interp(part.ener, spec_data.sig_ion_tot_data.first, spec_data.sig_ion_tot_data.second, true, false, 0., 0.);
  double sig_list[4] = {sig_scat, sig_brem, sig_exc, sig_ion};
  for ( size_t i = 0; i < 4; i++ ) {
    sig_tot += sig_list[i];
    sig_cum.push_back(sig_tot);
  }
  return findIdx(sig_tot*xi(), sig_cum);
}

/**
 * @brief Select an ionization.
 * 
 * @param Zelem The proton number of the element.
 * @return The ionization index.
*/
int Sim::choseIon(int Zelem) {
  Vector1d1dVector sig_ion_data = eedl[Zelem-1].sig_ion_data;
  double sig_tot = 0.0;
  Vector1d sig_cum;
  for ( size_t i = 0; i < sig_ion_data.size(); i++ ) {
    Vector1d1d sig_ion_data_1ion = sig_ion_data[i];
    double sig = interp(part.ener, sig_ion_data_1ion.first, sig_ion_data_1ion.second, true, false, 0., 0.);
    sig_tot += sig;
    sig_cum.push_back(sig_tot);
  }
  return findIdx(sig_tot*xi(), sig_cum);
}

/**
 * @brief Model an interaction.
 * 
 * Model an interacton:
 *  1. Select a species to interact with.
 *  2. If an element is selected, select an interaction type.
 *  3. If ionization is selected, select a particular ionization.
 * 
 * @param event The event object.
*/
void Sim::interact(Event &event) {
  event.Zelem = choseElem();
  if ( event.Zelem == 0 ) {
    event.interaction = flags::Bturb;
    part.newBvec(Bmag_co, Bmag_turb);
  } else {
    SpecData spec_data = eedl[event.Zelem-1];
    event.interaction = choseInter(event.Zelem);
    switch ( event.interaction ) {
      case flags::scat:
      event.cos_th = calcCosThScat(xi(), part.ener, spec_data.th_scat_data.first, spec_data.th_scat_data.second, spec_data.th_scat_data.third);
      part.scat(event.cos_th, part.vel);
      break;
      case flags::brem:
      event.ener_loss = interp(part.ener, spec_data.spec_brem_data.first, spec_data.spec_brem_data.second);
      part.loseEner(event.ener_loss);
      break;
      case flags::exc:
      event.ener_loss = interp(part.ener, spec_data.spec_exc_data.first, spec_data.spec_exc_data.second);
      part.loseEner(event.ener_loss);
      break;
      case flags::ion:
      event.ion = choseIon(event.Zelem);
      Vector1d2d2d spec_ion_data_1ion = spec_data.spec_ion_data[event.ion];
      event.ener_loss = calcEnerLoss(xi(), part.ener, spec_ion_data_1ion.first, spec_ion_data_1ion.second, spec_ion_data_1ion.third);
      event.ener_sec = event.ener_loss - spec_data.ener_bind_list[event.ion];
      part.loseEner(event.ener_loss);
      break;
    }
  }
  event.ener = part.ener;
}

/**
 * @brief Step the simulation forward one step.
*/
void Sim::step() {
  Event event(part.id, nstep);
  double sig_tot = calcSigTot();
  move(sig_tot, event);
  interact(event);
  event_list.push_back(event);
  nstep += 1;
}

/**
 * @brief Step the simulation forward multiple steps.
 * 
 * Step the simulation forward by a given number of steps.
 * Then write the event list to the outfile.
 * Finally, reset the event list.
 * 
 * @param num The number of steps.
*/
void Sim::multistep(int num) {
  for ( int i = 0; i < num; i++ ) { step(); }
  writeEvent(outfile, event_list);
  event_list.clear();
}

