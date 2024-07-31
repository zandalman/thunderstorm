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
Sim::Sim(Part part_, const EEDLData& eedl_, const Vector1d& ab_, std::string outfile_, double rho_, double temp_, double ion_state_avg_, double L_, double beta_, double mach_A_, double sig_turb_frac_, double alpha_, double cos_th_cut_)
  : part(part_)                   // The particle object.
  , eedl(eedl_)                   // Data from the EEDL database.
  , ab(ab_)                       // A vector of elemental abundances.
  , outfile(outfile_)             // The outfile name.
  , rho(rho_)                     // The density [g/cc].
  , temp(temp_)                   // The temperature [K].
  , ion_state_avg(ion_state_avg_) // The average ionization state.
  , L(L_)                         // The injection scale of the turbulence [cm].
  , beta(beta_)                   // The plasma beta
  , mach_A(mach_A_)               // The Alfven Mach number.
  , sig_turb_frac(sig_turb_frac_) // The effective turbulence cross section as a fraction of the total cross section.
  , alpha(alpha_)                 // The exponent of the magnetic field curvature spectrum
  , cos_th_cut(cos_th_cut_)       // The cutoff scattering angle cosine for discrete Moller scattering.
  , do_Bfield(beta > 0.)          // Whether a magnetic field is present.
  , do_turb(mach_A_ != 0.)        // Whether turbulence is present.
  , do_intermittancy(alpha_ > 0.) // Whether to include the effects of intermittancy in MHD turbulence.
  , nstep (0)                     // The step number.
  , time(0.0)                     // The simulation time [s].
  , n_i(0.0)                      // The ion number density [1/cc].
  , n_e_free(0.0)                 // The free electron number density [1/cc].
  , lam_deb(0.0)                  // The Debye length [1/cc].
  , B0(0.0)                       // The coherent magnetic field amplitude [G].
  , do_ion(ion_state_avg > 0.)    // Whether the atoms are ionized.
  , scale(0.)                     // The scale parameter for turbulent diffusion [cm^(-1/2)].
  , rperp_max(0.)                 // The truncation parameter for turbulent diffusion [cm].
  , Brms(0.)                      // The RMS magnetic field amplitude [G].

  {
    calcLamDeb(ab, rho, temp, ion_state_avg, n_i, n_e_free, lam_deb);
    calcStableParam(L, mach_A, scale, rperp_max);
    B0 = calcB0(n_i + n_e_free, temp, beta, mach_A);
    Brms = B0 * sqrt(1.0 + mach_A*mach_A);
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
  event_death.x = part.pos.x;
  event_death.y = part.pos.y;
  event_death.z = part.pos.z;
  event_death.interaction = flags::death;
  event_list.push_back(event_death);
  part.alive = false;
}

/**
 * @brief Compute the total cross section.
 * 
 * @return The total cross section [cm^2].
*/
double Sim::calcSigTot() {
  double sig_tot = 0.0;
  double sig_moller = do_ion ? calcSigMoller(part.gam(), part.beta(), lam_deb, cos_th_cut) : 0.;
  sig_tot += sig_moller * n_e_free / n_i;
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    double sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second, true, false, 0., 0.);
    sig_tot += sig * ab[i+1];
  }
  part.flag_intermittancy = false;
  double lam_intermittancy = (do_Bfield && do_turb && do_intermittancy) ? calcLamIntermittancy(part.m_i, part.q_i, part.gam(), part.vel.mag(), Brms, L, mach_A, alpha) : 0.0;
  double sig_intermittancy = (do_Bfield && do_turb && do_intermittancy) ? fabs(part.cos_alpha()) / (rho * constants::N_A * lam_intermittancy) : 0.0;
  if ( sig_intermittancy > sig_tot ) {
    part.flag_intermittancy = true;
    part.lam_intermittancy = lam_intermittancy;
    sig_intermittancy = 0.0;
  }
  sig_tot += sig_intermittancy;
  double sig_turb = (do_Bfield && do_turb) ? sig_turb_frac * sig_tot / (1.0 - sig_turb_frac) : 0.;
  sig_tot += sig_turb;
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
    if ( part.flag_intermittancy ) {
      dis = calcIntermittancyTrans(xi(), xi(), part.lam_intermittancy, dis, part.cos_alpha());
    }
    part.pos = part.pos + dis * part.cos_alpha() * part.Bvec.unit();
    if ( part.flag_turb_diff ) {
      double lam_turb = fabs(part.cos_alpha()) / (rho * constants::N_A * sig_tot * sig_turb_frac);
      double rperp = calcTurbDiff(xi(), xi(), lam_turb, scale, rperp_max);
      part.pos = part.pos + randPerpVec(part.Bvec, rperp);
      part.flag_turb_diff = false;
    }
    part.vel = rotate(part.vel, part.Bvec, cos(2.*M_PI * xi()));
    event.ener_loss_sync = calcPowerSync(part.m_i, part.q_i, part.gam(), part.beta(), part.Bvec.mag(), part.cos_alpha()) * dt;
  } else {
    part.pos = part.pos + dis * part.vel.unit();
  }
  if ( do_ion ) {
    event.ener_loss_cher = calcPowerCher(part.beta(), temp, n_e_free) * dt;
    event.ener_loss_moller = calcPowerMoller(part.ener, part.gam(), part.beta(), n_e_free, lam_deb, cos_th_cut) * dt;
  }
  part.loseEner(event.ener_loss_sync + event.ener_loss_cher + event.ener_loss_moller);
  event.time = time;
  event.x = part.pos.x;
  event.y = part.pos.y;
  event.z = part.pos.z;
}

/**
 * @brief Select an element.
 * 
 * @return The proton number of the selected element, or a flag for a non-element interactions.
*/
int Sim::choseElem() {
  if ( do_Bfield && do_turb && ( xi() < sig_turb_frac ) ) {
    return flags_elem::turb;
  }
  double sig_tot = 0.0;
  Vector1d sig_cum;
  double sig_moller = do_ion ? calcSigMoller(part.gam(), part.beta(), lam_deb, cos_th_cut) : 0.;
  sig_tot += sig_moller * n_e_free / n_i;
  sig_cum.push_back(sig_tot);
  double sig_intermittancy = (do_Bfield && do_turb && do_intermittancy && !part.flag_intermittancy) ? fabs(part.cos_alpha()) / (rho * constants::N_A * part.lam_intermittancy) : 0.0;
  sig_tot += sig_intermittancy;
  sig_cum.push_back(sig_tot);
  for ( size_t i = 0; i < eedl.size(); i++ ) {
    SpecData spec_data = eedl[i];
    double sig = interp(part.ener, spec_data.sig_tot_data.first, spec_data.sig_tot_data.second, true, false, 0., 0.);
    sig_tot += sig * ab[i+1];
    sig_cum.push_back(sig_tot);
  }
  int idx_elem = findIdx(sig_tot * xi(), sig_cum);
  switch ( idx_elem ) { 
    case 0:
    return flags_elem::moller;
    case 1:
    return flags_elem::intermittancy;
    default:
    return idx_elem - 1;
  }
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
  return findIdx(sig_tot*xi(), sig_cum) + 1;
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
  event.Zelem = choseElem();
  switch ( event.Zelem ) {
    case flags_elem::turb:
    event.interaction = flags::turb;
    part.flag_turb_diff = true;
    break;
    case flags_elem::moller:
    event.interaction = flags::moller;
    calcCosThScatEnerLossMoller(xi(), part.ener, part.gam(), part.beta(), lam_deb, cos_th_cut, event.cos_th, event.ener_loss);
    part.scat(event.cos_th, part.vel);
    part.loseEner(event.ener_loss);
    break;
    case flags_elem::intermittancy:
    event.interaction = flags::intermittancy;
    part.vel = randVec(part.vel.mag());
    break;
    default:
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
      event.ener_sec = calcEnerLoss(xi(), part.ener, spec_ion_data_1ion.first, spec_ion_data_1ion.second, spec_ion_data_1ion.third);
      event.ener_loss = event.ener_sec + spec_data.ener_bind_list[event.ion];
      part.loseEner(event.ener_loss);
      break;
    }
    break;
  }
  event.ener = part.ener;
  event.cos_alpha = part.cos_alpha();
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
