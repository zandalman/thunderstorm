// includes
#include <cmath>
#include <algorithm>
#include <iostream>

// headers
#include "const.h"
#include "random.h"
#include "part.h"
#include "sim.h"
#include "parser.h"
#include "functions.h"
#include "io.h"

/// @brief A constructor to initial to the Sim structure.
Sim::Sim(
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
)
  : part(part_)                   // The particle object.
  , eedl(eedl_)                   // Data from the EEDL database.
  , ab(ab_)                       // A vector of elemental abundances.
  , outfile(outfile_)             // The outfile name.
  , rho(rho_)                     // The density [g/cc].
  , temp(temp_)                   // The temperature [K].
  , ion_state(ion_state_)         // The ionization state vector.
  , q_avg(0.0)                    // The average ionization state.
  , qsq_avg(0.0)                  // The average square ionization state.
  , cos_th_cut(cos_th_cut_)       // The cutoff scattering angle cosine for discrete Moller scattering.
  , nstep (0)                     // The step number.
  , time(0.0)                     // The simulation time [s].
  , n_i(0.0)                      // The ion number density [1/cc].
  , n_e_free(0.0)                 // The free electron number density [1/cc].
  , lam_deb(0.0)                  // The Debye length [1/cc].
  , B0(B0_)                       // The coherent magnetic field amplitude [G].
  , neutral(neutral_)             // Whether the ejecta is neutral.
  , do_cerenkov(do_cerenkov_)     // Do Cerenkov energy losses.
  , do_sync(do_sync_)             // Do synchrotron energy losses.
  , mmw(mmw_)                     // Mean molecular weight [g/mol].
  {
    normIonState(ion_state, q_avg, qsq_avg);
    if ( !neutral ) calcLamDeb(ab, rho, temp, q_avg, qsq_avg, n_i, n_e_free, lam_deb);
  }

/**
 * @brief Reset the simulation with a new particle.
 * 
 * @param new_part The new particle.
 */
void Sim::reset(Part new_part) {
  part = new_part;
  nstep = 0;
  time = 0.0;
}

/**
 * @brief Kill the particle.
 */
void Sim::kill() { 
  Event event_death(part.id, nstep);
  event_death.time = time;
  event_death.splus = part.splus;
  event_death.sminus = part.sminus;
  event_death.ener = part.ener;
  event_death.cos_alpha = part.cos_alpha;
  event_death.interaction = flags::death;
  event_list.push_back(event_death);
  part.alive = false;
}

/**
 * @brief Compute the total cross section.
 * 
 * @return The total cross section [cm^2 mol / g].
*/
double Sim::calcSigTot() {
  double sig;
  double sig_tot = 0.0;
  double sig_moller = !neutral ? calcSigMoller(part.gam(), part.beta(), lam_deb, cos_th_cut) : 0.;
  sig_tot += sig_moller * n_e_free / n_i / mmw;
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    for ( size_t j = 0; j < ion_state.size(); j++ ) {
      if ( ion_state[j] == 0.0 ) continue;
      sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second[j], true, false, 0., 0.);
      sig_tot += sig * ab[i+1] * ion_state[j];
    }
  }
  return sig_tot;
}

/**
 * @brief Transport the particle.
 * 
 * @param sig_tot The total cross section [cm^2 mol / g].
 * @param event   The event structure.
*/
void Sim::move(double sig_tot, Event &event) {
  // compute the distance and travel time
  double dis = -log(1 - xi()) / (rho * constants::N_A * sig_tot);
  double dt = dis / (constants::c * part.beta());
  time += dt;
  // move the particle along the field line
  if ( part.cos_alpha >= 0 ) {
    part.splus += dis * part.cos_alpha;
  } else {
    part.sminus += -dis * part.cos_alpha;
  }
  // calculate energy loss and pitch angle diffusion in transport
  double edot_moller, cos_al_dot_sq_moller, cos_al_dot_sq_mott;
  if ( do_sync ) event.ener_loss_sync = calcPowerSync(part.m_i, part.q_i, part.gam(), part.beta(), B0, part.cos_alpha) * dt;
  if ( do_cerenkov && !neutral ) event.ener_loss_cher = calcPowerCerenkov(part.beta(), temp, n_e_free) * dt;
  if ( !neutral ) {
    calcPowerDiffMoller(
      part.ener, 
      part.gam(), 
      part.beta(), 
      part.cos_alpha, 
      n_e_free, 
      lam_deb, 
      cos_th_cut, 
      edot_moller, 
      cos_al_dot_sq_moller
    );
    calcDiffMott(
      part.ener,
      part.gam(),
      part.beta(),
      part.cos_alpha,
      n_i,
      lam_deb,
      qsq_avg,
      mmw,
      eedl,
      ab,
      cos_al_dot_sq_mott
    );
    event.ener_loss_moller = edot_moller * dt;
    event.dlt_cos_al_moller = sqrt(cos_al_dot_sq_moller) * dt;
    event.dlt_cos_al_mott = sqrt(cos_al_dot_sq_mott) * dt;
  }

  part.loseEner(event.ener_loss_sync + event.ener_loss_cher + event.ener_loss_moller);
  part.scatCum(sqrt(cos_al_dot_sq_moller + cos_al_dot_sq_mott) * dt);

  // update event
  event.time = time;
  event.splus = part.splus;
  event.sminus = part.sminus;
}

/**
 * @brief Select an element.
 * 
 * @param moller Whether the interaction is Moller scattering.
 * @param Zelem  The proton number of the selected element.
 * @param stage  The selected ionization stage.
*/
void Sim::choseElem(bool &moller, int &Zelem, int &stage) {
  double sig;
  double sig_tot = 0.0;
  Vector1d sig_tot1d;
  double sig_tot_1elem;
  Vector1d sig_cum1d;
  Vector2d sig_cum2d;
  Vector1d sig_cum2d_1elem;
  double sig_moller = neutral ? 0.0: calcSigMoller(part.gam(), part.beta(), lam_deb, cos_th_cut);
  sig_tot += sig_moller * n_e_free / n_i / mmw;
  sig_cum1d.push_back(sig_tot);
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    sig_tot_1elem = 0.0;
    sig_cum2d_1elem.clear();
    for ( size_t j = 0; j < ion_state.size(); j++ ) {
      sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second[j], true, false, 0., 0.);
      sig_tot += sig * ab[i+1] * ion_state[j];
      sig_tot_1elem += sig * ab[i+1] * ion_state[j];
      sig_cum2d_1elem.push_back(sig_tot_1elem);
    }
    sig_tot1d.push_back(sig_tot_1elem);
    sig_cum1d.push_back(sig_tot);
    sig_cum2d.push_back(sig_cum2d_1elem);
  }
  int idx_elem = findIdx(sig_tot * xi(), sig_cum1d);
  if ( idx_elem == 0 ) {
    moller = true;
    Zelem = -1;
    stage = -1;
  } else {
    moller = false;
    Zelem = idx_elem;
    stage = findIdx(sig_tot1d[Zelem - 1] * xi(), sig_cum2d[Zelem - 1]);
  }
}

/**
 * @brief Select an interaction.
 * 
 * @param Zelem The proton number of the element.
 * @param stage The ionization stage of the element.
 * @return The interaction flag.
*/
int Sim::choseInter(int Zelem, int stage) {
  SpecData spec_data = eedl[Zelem-1];
  double sig_tot = 0.0;
  Vector1d sig_cum;
  double sig_scat = interp(part.ener, spec_data.sig_scat_data.first, spec_data.sig_scat_data.second, true, false, 0., 0.);
  double sig_brem = interp(part.ener, spec_data.sig_brem_data.first, spec_data.sig_brem_data.second, true, false, 0., 0.);
  double sig_exc = interp(part.ener, spec_data.sig_exc_data.first, spec_data.sig_exc_data.second, true, false, 0., 0.);
  double sig_ion = interp(part.ener, spec_data.sig_ion_tot_data.first, spec_data.sig_ion_tot_data.second[stage], true, false, 0., 0.);
  double sig_list[4] = {sig_scat, sig_brem, sig_exc, sig_ion};
  for ( size_t i = 0; i < 4; i++ ) {
    sig_tot += sig_list[i];
    sig_cum.push_back(sig_tot);
  }
  return findIdx(sig_tot * xi(), sig_cum) + 1;
}

/**
 * @brief Select an ionization.
 * 
 * @param Zelem The proton number of the element.
 * @param stage The ionization stage of the element.
 * @return The ionization index.
*/
int Sim::choseIon(int Zelem, int stage) {
  Vector1d1dVector1d sig_ion_data = eedl[Zelem-1].sig_ion_data[stage];
  double sig_tot = 0.0;
  Vector1d sig_cum;
  for ( size_t i = 0; i < sig_ion_data.size(); i++ ) {
    Vector1d1d sig_ion_data_1ion = sig_ion_data[i];
    double sig = interp(part.ener, sig_ion_data_1ion.first, sig_ion_data_1ion.second, true, false, 0., 0.);
    sig_tot += sig;
    sig_cum.push_back(sig_tot);
  }
  return findIdx(sig_tot * xi(), sig_cum);
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
  bool moller;
  int idx_ion;
  choseElem(moller, event.Zelem, event.stage);
  if ( moller ) {
    event.interaction = flags::moller;
    calcCosThScatEnerLossMoller(xi(), part.ener, part.gam(), part.beta(), lam_deb, cos_th_cut, event.cos_th, event.ener_loss);
    part.scat(xi(), event.cos_th);
    part.loseEner(event.ener_loss);
  } else {
    SpecData spec_data = eedl[event.Zelem - 1];
    event.interaction = choseInter(event.Zelem, event.stage);
    switch ( event.interaction ) {
      case flags::scat:
      event.cos_th = calcCosThScat(xi(), part.ener, spec_data.th_scat_data.first, spec_data.th_scat_data.second, spec_data.th_scat_data.third);
      part.scat(xi(), event.cos_th);
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
      idx_ion = choseIon(event.Zelem, event.stage);
      event.ion = spec_data.ss_list[idx_ion];
      Vector1d2d2d spec_ion_data_1ion = spec_data.spec_ion_data[event.stage][idx_ion];
      double ener_bind = spec_data.ener_bind_list[event.stage][idx_ion];
      event.ener_sec = calcEnerLoss(xi(), part.ener, part.ener - ener_bind, spec_ion_data_1ion.first, spec_ion_data_1ion.second, spec_ion_data_1ion.third);
      event.ener_loss = event.ener_sec + ener_bind;
      part.loseEner(event.ener_loss);
      break;
    }
  }
  event.ener = part.ener;
  event.cos_alpha = part.cos_alpha;
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
