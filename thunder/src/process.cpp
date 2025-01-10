// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <cstdio>
#include <mpi.h>

// headers
#include "const.h"
#include "functions.h"
#include "parser.h"
#include "process.h"
#include "io.h"
#include "random.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;
template <typename T>
using vector3d = std::vector<vector2d<T>>;

/// @brief A constructor to initialize the Data structure.
Data::Data(
  double mach_A_, 
  double scale_, 
  double ener_, 
  MiscParam misc_param_,
  const std::vector<Stat> &stat_list
  )
  : mach_A(mach_A_)
  , ener(ener_)
  , ener_min(misc_param_.ener_min)
  , scale(scale_)
  , turb(misc_param_.turb * scale_)
  , escaped(false)
  , ener_start(ener_)
  , time_start(0.0)
  , ener_prev(ener_)
  , time_prev(0.0)
  , splus_prev(0.0)
  , sminus_prev(0.0)
  , lam_scat(0.0)
  , s_scat(0.0)
  , oss()
 {
  Stat stat;
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    stat = stat_list[i];
    part_stat_list.push_back(std::vector<double>(stat.size, 0.0));
    mean_stat_list.push_back(std::vector<double>(stat.size, 0.0));
    M2_stat_list.push_back(std::vector<double>(stat.size, 0.0));
    M3_stat_list.push_back(std::vector<double>(stat.size, 0.0));
    M4_stat_list.push_back(std::vector<double>(stat.size, 0.0));
  }
  lam_scat = mach_A > 1.0 ? turb / (mach_A*mach_A*mach_A) : turb * mach_A*mach_A*mach_A*mach_A;
  s_scat = -log(1.0 - xi()) * lam_scat;
  pos = Vec(0.0, 0.0, 0.0);
  Bhat = calcRandVec(mach_A);
 }

/// @brief Reset the particle data.
void Data::reset() {
  escaped = false;
  ener_start = ener;
  time_start = 0.0;
  ener_prev = ener;
  time_prev = 0.0;
  splus_prev = 0.0;
  sminus_prev = 0.0;
  s_scat = -log(1.0 - xi()) * lam_scat;
  pos = Vec(0.0, 0.0, 0.0);
  Bhat = calcRandVec(mach_A);
  oss.str(""); oss.clear();
  
  for ( size_t i = 0; i < part_stat_list.size(); i++ ) {
    std::fill(part_stat_list[i].begin(), part_stat_list[i].end(), 0.0);
  }
}

/**
 * @brief Update the statistics with a new particle.
 * 
 * @param n_int     The current particle count.
 * @param stat_list The list of statistics.
 */
void Data::calcStat(int n_int, const std::vector<Stat> &stat_list) {
  double delta, delta_np1, delta_np1_sq, term1;
  double n = static_cast<double>(n_int);
  double np1 = n + 1.0;
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    const Stat &stat = stat_list[i];
    for ( size_t j = 0; j < stat.size; j++ ) {
      if ( n_int == 0 ) {
        mean_stat_list[i][j] = part_stat_list[i][j];
      } else {
        delta = part_stat_list[i][j] - mean_stat_list[i][j];
        delta_np1 = delta / np1;
        delta_np1_sq = delta_np1*delta_np1;
        term1 = delta * delta_np1 * n;
        M4_stat_list[i][j] += term1 * delta_np1_sq * (np1*np1 - 3.0 * np1 + 3.0) + 6.0 * delta_np1_sq * M2_stat_list[i][j] - 4.0 * delta_np1 * M3_stat_list[i][j];
        M3_stat_list[i][j] += term1 * delta_np1 * (np1 - 2.0) - 3.0 * delta_np1 * M2_stat_list[i][j];
        M2_stat_list[i][j] += term1;
        mean_stat_list[i][j] += delta_np1;
      }
    }
  }
}

/**
 * @brief Process a binary data file output by thunderstorm.
 * 
 * @param datafile_name       The name of the binary data file.
 * @param num_event_per_chunk The number of events per chunk.
 * @param bin_list            The list of bins.
 * @param stat_list           The list of statistics.
 * @param count               The particle count.
 * @param data_grid           The grid of data.  
 */
void processFile(
  const std::string &datafile_name, 
  size_t num_event_per_chunk, 
  const vector2d<double> &bin_list,
  const std::vector<Stat>& stat_list, 
  const std::string &histdir_name,
  int idx_hist_max,
  int &idx_hist,
  int &count, 
  vector3d<Data>& data_grid
) {

  // Open file
  std::ifstream datafile(datafile_name, std::ios::binary);
  if ( !datafile.is_open() ) {
    std::cerr << "Failed to open file " << datafile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Compute chunk size
  const size_t event_size = sizeof(Event);
  const size_t chunk_size = event_size * num_event_per_chunk;

  // Allocate memory to hold the chunk content
  std::vector<char> buffer(chunk_size);
  Event* event;
  
  // Read data in chunks
  int current_id = -1;
  std::ostringstream oss;
  bool do_hist;

  while ( datafile.read(buffer.data(), chunk_size) ) {
    
    for ( size_t i = 0; i < num_event_per_chunk; i++ ) {
      do_hist = idx_hist < idx_hist_max;
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if ( current_id == -1 ) { current_id = event->id; }
      if ( current_id == event->id ) {
        processEvent(event, do_hist, bin_list, data_grid);      
      } else {
        if ( do_hist ) {
          oss << histdir_name << "/hist";
          oss << std::setw(5) << std::setfill('0') << idx_hist << ".txt";
          clearFile(oss.str());
          writeHist(oss.str(), bin_list, data_grid);
          oss.str(""); oss.clear();
          idx_hist++;
        }
        postProcPart(count, stat_list, data_grid);
        current_id = event->id;
        count++;
      }
    }
  }

  // Handle the case where the last chunk might not be full
  if ( datafile.eof() ) {
    size_t bytes_read = datafile.gcount();
    if ( bytes_read > 0 ) {
      size_t num_event_last_chunk = bytes_read / event_size;
      for ( size_t i = 0; i < num_event_last_chunk; i++ ) {
        do_hist = idx_hist < idx_hist_max;
        event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
        if ( current_id == -1 ) { current_id = event->id; }
        if ( current_id == event->id ) {
          processEvent(event, do_hist, bin_list, data_grid);
        } else {
          if ( do_hist ) {
            oss << histdir_name << "/hist";
            oss << std::setw(5) << std::setfill('0') << idx_hist << ".txt";
            clearFile(oss.str());
            writeHist(oss.str(), bin_list, data_grid);
            oss.str(""); oss.clear();
            idx_hist++;
          }
          postProcPart(count, stat_list, data_grid);
          current_id = event->id;
          count++;
        }
      }
    }
  } else {
    std::cerr << "Error reading file " << datafile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  datafile.close();
}

/**
 * @brief Post-process event data for a single particle and reset the data grid.
 * 
 * @param count     The particle count.
 * @param stat_list The list of statistics.
 * @param data_grid The grid of data.
 */
void postProcPart(
  int count,
  const std::vector<Stat> &stat_list, 
  vector3d<Data> &data_grid
) {
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        Data &data = data_grid[i][j][k];
        
        // compute thermalization efficiency
        if ( data.escaped ) {
          data.part_stat_list[stat_tag::eps_thm][0] = fmax(0.0, 1.0 - data.ener_prev / data.ener);
        } else {
          data.part_stat_list[stat_tag::eps_thm][0] = 1.0;
        }
        
        // aggregate statistics and reset
        data.calcStat(count, stat_list);
        data.reset();
      }
    }
  }
}

/**
 * @brief Process an event and add the post-processed data to the particle data.
 * 
 * @param event     The event struct.
 * @param bin_list  The list of bins.
 * @param data_grid The grid of data.
 */
void processEvent(
  const Event* event, 
  bool do_hist,
  const vector2d<double> &bin_list, 
  vector3d<Data>& data_grid
) {
  // throw error if event pointer is null.
  if ( event == nullptr ) {
    std::cerr << "Error: Null event pointer." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int flag;
  double ener_loss, time_rel;
  double dt, dsplus, dsminus, ds, sign;
  size_t idx_ener_sec, idx_time, idx_ener;
  
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        
        Data &data = data_grid[i][j][k];
        if ( data.escaped ) continue;
        
        // if energy is above starting energy, update start time and coordinates
        if ( event->ener > data.ener ) {
          data.ener_start = event->ener;
          data.time_start = event->time;
          data.ener_prev = event->ener;
          data.time_prev = event->time;
          data.splus_prev = event->splus;
          data.sminus_prev = event->sminus;
          continue;
        }

        // compute useful quantities
        ener_loss = data.ener_prev - event->ener;
        time_rel = event->time - data.time_start;
        idx_ener = findIdx(event->ener, bin_list[bin_tag::ener]);
        idx_time = findIdx(time_rel, bin_list[bin_tag::time]);

        // compute transport
        dt = event->time - data.time_prev;
        dsplus = event->splus - data.splus_prev;
        dsminus = event->sminus - data.sminus_prev;
        ds = dsplus + dsminus;
        sign = dsplus > dsminus ? 1.0 : -1.0;

        while ( true ) {
          if ( ds < data.s_scat ) {
            data.pos = data.pos + sign * ds * data.Bhat;
            data.s_scat = data.s_scat - ds;
            break;
          } else {
            data.pos = data.pos + sign * data.s_scat * data.Bhat;
            ds = ds - data.s_scat;
            data.s_scat = -log(1.0 - xi()) * data.lam_scat;
            data.Bhat = calcRandVec(data.mach_A); // Resample B-field direction
          }
        }

        // update time and distance
        data.ener_prev = event->ener;
        data.time_prev = event->time;
        data.splus_prev = event->splus;
        data.sminus_prev = event->sminus;

        // compute energy histograms
        if ( idx_ener > 0 && idx_ener < bin_list[bin_tag::ener].size() ) {
          data.part_stat_list[stat_tag::time_ener][idx_ener - 1] += dt;
        }
        
        // compute time histograms
        if ( idx_time > 0 && idx_time < bin_list[bin_tag::time].size() ) {
          data.part_stat_list[stat_tag::ener_loss_time][idx_time - 1] += ener_loss;
        }
        
        // compute interaction histograms
        flag = event->interaction;
        switch ( flag ) {
          case flags::scat: // scattering
          data.part_stat_list[stat_tag::num_ev_inter][inter_tag::scat] += 1.0;
          break;
          case flags::brem: // Bremsstrahlung
          data.part_stat_list[stat_tag::num_ev_inter][inter_tag::brem] += 1.0;
          data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::brem] += event->ener_loss;
          break;
          case flags::exc: // excitation
          data.part_stat_list[stat_tag::num_ev_inter][inter_tag::exc] += 1.0;
          data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::exc] += event->ener_loss;
          break;
          case flags::ion: // ionization
          data.part_stat_list[stat_tag::num_ev_inter][inter_tag::ion] += 1.0;
          data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::ion] += event->ener_loss;
          data.part_stat_list[stat_tag::num_ion_elem][event->Zelem - 1] += 1.0;
          // compute secondary energy histograms
          idx_ener_sec = findIdx(event->ener_sec, bin_list[bin_tag::ener_sec]);
          if (idx_ener_sec > 0 && idx_ener_sec < bin_list[bin_tag::ener_sec].size()) {
            data.part_stat_list[stat_tag::num_sec_ener][idx_ener_sec - 1] += 1.0;
          }
          break;
          case flags::moller: // Moller
          data.part_stat_list[stat_tag::num_ev_inter][inter_tag::moller] += 1.0;
          if ( !std::isnan(event->ener_loss) ) { // for some reason, this is sometimes NaN
            data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::moller] += event->ener_loss;
          }
          break;
        }
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::moller] += event->ener_loss_moller;
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::sync] += event->ener_loss_sync;
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::cher] += event->ener_loss_cher;

        // check if particle has crossed thermalization or escape barrier
        if ( fabs(data.pos.z) > 0.5 * data.scale ) {
          data.escaped = true;
          flag = -1;
          if (idx_ener > 0 && idx_ener < bin_list[bin_tag::ener].size()) {
            data.part_stat_list[stat_tag::num_escape][idx_ener - 1] += 0.5;
          }
        }

        // write data
        if ( do_hist ) {
          data.oss << time_rel << ",";
          data.oss << data.pos.x << "," << data.pos.y << "," << data.pos.z << ",";
          data.oss << event->cos_alpha << ", ";
          data.oss << event->ener << ", ";
          data.oss << flag << std::endl;
        }
      }
    }
  }
}

/**
 * @brief Flatten the data so it can be communicated via MPI.
 * 
 * @param data_grid           The grid of data.
 * @param stat_list           The list of statistics.
 * @param mean_stat_list_flat The flattened list of mean statistics.
 * @param M2_stat_list_flat   The flattened list of M2 statistics.
 * @param M3_stat_list_flat   The flattened list of M3 statistics.
 * @param M4_stat_list_flat   The flattened list of M4 statistics.
 * @param size_t              The size of the flattened list of statistics.
 */
void getFlatData(
  const vector3d<Data>& data_grid,
  const std::vector<Stat> &stat_list, 
  std::vector<double> &mean_stat_list_flat, 
  std::vector<double> &M2_stat_list_flat, 
  std::vector<double> &M3_stat_list_flat, 
  std::vector<double> &M4_stat_list_flat, 
  size_t &size_flat
) {
  size_flat = 0;
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        const Data &data = data_grid[i][j][k];
        for ( size_t l = 0; l < stat_list.size(); l++ ) {
          const Stat &stat = stat_list[l];
          for ( size_t m = 0; m < stat.size; m++ ) {
            mean_stat_list_flat.push_back(data.mean_stat_list[l][m]);
            M2_stat_list_flat.push_back(data.M2_stat_list[l][m]);
            M3_stat_list_flat.push_back(data.M2_stat_list[l][m]);
            M4_stat_list_flat.push_back(data.M2_stat_list[l][m]);
            size_flat += 1;
          }
        }
      }
    }
  }
}
