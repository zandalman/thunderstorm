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
#include <chrono>
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
  double dx_, 
  double ener_low_, 
  double ener_high_,
  double ener_min_,
  double dt_,
  int ndim_,
  int nmom_,
  const std::vector<Stat> &stat_list
  )
  
  : mach_A(mach_A_)
  , dx(dx_)
  , ener_low(ener_low_)
  , ener_high(ener_high_)
  , ener_min(ener_min_)
  , dt(dt_)
  , super(mach_A_ >= 1.0)
  , ndim(ndim_)
  , nker(2 * ndim_)
  , nmom(nmom_)
  , nstat(stat_list.size())

  , ener(0.0)
  , ener_prev(0.0)
  , ener_start(0.0)

  , time(0.0)
  , time_prev(0.0)
  , time_start(0.0)

  , outoftime(false)
  , thermalized(false)

  , splus_prev(0.0)
  , sminus_prev(0.0)
  , lam_scat(0.0)
  , s_scat(0.0)
  , oss()
 {
  size_t size_stat;
  part_stat_list.resize(nstat);
  M1_stat_list.resize(nstat);
  if (nmom >= 2) M2_stat_list.resize(nstat);
  if (nmom >= 3) M3_stat_list.resize(nstat);
  if (nmom >= 4) M4_stat_list.resize(nstat);
  for ( size_t i = 0; i < nstat; i++ ) {
    size_stat = stat_list[i].size;
    part_stat_list[i].resize(nker);
    M1_stat_list[i].resize(nker);
    if (nmom >= 2) M2_stat_list[i].resize(nker);
    if (nmom >= 3) M3_stat_list[i].resize(nker);
    if (nmom >= 4) M4_stat_list[i].resize(nker);
    for ( size_t j = 0; j < nker; j++ ) {
      part_stat_list[i][j].resize(size_stat, 0.0);
      M1_stat_list[i][j].resize(size_stat, 0.0);
      if (nmom >= 2) M2_stat_list[i][j].resize(size_stat, 0.0);
      if (nmom >= 3) M3_stat_list[i][j].resize(size_stat, 0.0);
      if (nmom >= 4) M4_stat_list[i][j].resize(size_stat, 0.0);
    }
  }
  lam_scat = super ? dx * mach_A*mach_A*mach_A : dx;
  
  // initialize the energy to a random energy within the bin
  // initial energies are log-spaced within each bin to match the bin spacing
  ener = ener_low * pow(ener_high / ener_low, xi());
  ener_start = ener;
  ener_prev = ener;

  // initialize the time to a random time within the timestep
  time = xi() * dt;
  time_start = time;
  time_prev = time;
  outoftime = false;
  thermalized = false;

  splus_prev = 0.0;
  sminus_prev = 0.0;

  s_scat = -log(1.0 - xi()) * lam_scat;
  pos = Vec(xi(), xi(), xi()) * dx;
  Bhat = calcRandVec(mach_A, super);

  oss.str(""); oss.clear();

 }

/// @brief Reset the particle data.
void Data::reset() {
  
  // initialize the energy to a random energy within the bin
  // initial energies are log-spaced within each bin to match the bin spacing
  ener = ener_low * pow(ener_high / ener_low, xi());
  ener_start = ener;
  ener_prev = ener;

  // initialize the time to a random time within the timestep
  time = xi() * dt;
  time_start = time;
  time_prev = time;
  outoftime = false;
  thermalized = false;

  splus_prev = 0.0;
  sminus_prev = 0.0;

  s_scat = -log(1.0 - xi()) * lam_scat;
  pos = Vec(xi(), xi(), xi()) * dx;
  Bhat = calcRandVec(mach_A, super);

  oss.str(""); oss.clear();
  
  for ( size_t i = 0; i < nstat; i++ ) {
    for ( size_t j = 0; j < nker; j++ ) {
      std::fill(part_stat_list[i][j].begin(), part_stat_list[i][j].end(), 0.0);
    }
  }
}

/**
 * @brief Update the statistics with a new particle.
 * See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
 * 
 * @param n_int     The current particle count.
 * @param stat_list The list of statistics.
 */
void Data::calcStat(int n_int, const std::vector<Stat> &stat_list) {
  double delta, delta_np1, delta_np1_sq, term1;
  double part_stat, M1, M2, M3;
  double n = static_cast<double>(n_int);
  double np1 = n + 1.0;
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    const Stat &stat = stat_list[i];
    for ( size_t j = 0; j < nker; j++ ) {
      for ( size_t k = 0; k < stat.size; k++ ) {
        part_stat = part_stat_list[i][j][k];
        M1 = M1_stat_list[i][j][k];
        if ( nmom >= 2 ) M2 = M2_stat_list[i][j][k];
        if ( nmom >= 3 ) M3 = M3_stat_list[i][j][k];
        if ( n_int == 0 ) {
          M1_stat_list[i][j][k] = part_stat;
        } else {
          delta = part_stat - M1;
          delta_np1 = delta / np1;
          delta_np1_sq = delta_np1*delta_np1;
          term1 = delta * delta_np1 * n;
          if ( nmom >= 4 ) M4_stat_list[i][j][k] += term1 * delta_np1_sq * (np1*np1 - 3.0 * np1 + 3.0) + 6.0 * delta_np1_sq * M2 - 4.0 * delta_np1 * M3;
          if ( nmom >= 3 ) M3_stat_list[i][j][k] += term1 * delta_np1 * (np1 - 2.0) - 3.0 * delta_np1 * M2;
          if ( nmom >= 2 ) M2_stat_list[i][j][k] += term1;
          M1_stat_list[i][j][k] += delta_np1;
        }
      }
    }
  }
}

/**
 * @brief Process a binary data file output by thunderstorm.
 * 
 * @param datafile_name       The name of the binary data file.
 * @param num_event_per_chunk The number of events per chunk.
 * @param start               The start time of the simulation.
 * @param walltime            The wallclock time of the job [s].
 * @param bin_list            The list of bins.
 * @param stat_list           The list of statistics.
 * @param count               The particle count.
 * @param data_grid           The grid of data.  
 * @param no_time             Whether the job has run out of wallclock time.
 */
void processFile(
  const std::string &datafile_name, 
  size_t num_event_per_chunk, 
  std::chrono::steady_clock::time_point start,
  int walltime,
  const vector2d<double> &bin_list,
  const std::vector<Stat>& stat_list, 
  const std::string &histdir_name,
  int idx_hist_max,
  int &idx_hist,
  int &count, 
  vector3d<Data>& data_grid,
  bool &no_time
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

  // Define lambda to process chunks
  auto processChunk = [&](size_t num_event) {
    for ( size_t i = 0; i < num_event; i++ ) {
      do_hist = idx_hist < idx_hist_max;
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if (current_id == -1) current_id = event->id;
      if (current_id == event->id) {
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
        postProcPart(count, bin_list, stat_list, data_grid);
        current_id = event->id;
        count++;
      }
    }
  };

  while ( datafile.read(buffer.data(), chunk_size) ) {
    
    // start timer
    auto start_chunk = std::chrono::steady_clock::now();
    
    // process chunk
    processChunk(num_event_per_chunk);

    // break if we are about to exceed walltime
    auto now = std::chrono::steady_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    auto chunktime = std::chrono::duration_cast<std::chrono::seconds>(now - start_chunk).count();
    if ( runtime > walltime - 2.0 * chunktime ) {
      no_time = true;
      break;
    }
  }

  // Handle the case where the last chunk might not be full
  if ( !no_time && datafile.eof() ) {
    size_t bytes_read = datafile.gcount();
    if ( bytes_read > 0 ) {
      
      // compute number of events in last chunk
      size_t num_event_last_chunk = bytes_read / event_size;
      
      // process last chunk
      processChunk(num_event_last_chunk);
    }
  } else if ( !no_time ) {
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
  const vector2d<double> &bin_list,
  const std::vector<Stat> &stat_list, 
  vector3d<Data> &data_grid
) {
  
  size_t iker;
  int idx_ener;
  
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        
        Data &data = data_grid[i][j][k];
        iker = calcKer(data.pos, data.dx, data.ndim);
        
        if ( data.outoftime ) {
          idx_ener = findIdx(data.ener_prev, bin_list[bin_tag::ener]);
          if ( idx_ener > 0 && idx_ener < bin_list[bin_tag::ener].size() ) {
            data.part_stat_list[stat_tag::ener][iker][idx_ener - 1] = data.ener_prev / data.ener_start;
          }
        } else if ( data.thermalized ) {
          data.part_stat_list[stat_tag::ener_thm][iker][0] += data.ener_prev / data.ener_start;
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

  int flag, iker;
  double ener_loss, time;
  double dt, dsplus, dsminus, ds, sign;
  size_t idx_ener_sec, idx_ener;
  
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        
        Data &data = data_grid[i][j][k];
        if ( data.outoftime || data.thermalized ) continue;
        
        // if energy is above starting energy, update initial time and coordinates
        if ( event->ener > data.ener ) {
          data.ener_start = event->ener;
          data.time_start = event->time;
          data.ener_prev = event->ener;
          data.time_prev = event->time;
          data.splus_prev = event->splus;
          data.sminus_prev = event->sminus;
          continue;
        }

        // compute relative time and energy loss
        time = data.time + event->time - data.time_start;
        ener_loss = data.ener_prev - event->ener;

        // compute kernel flag
        iker = calcKer(data.pos, data.dx, data.ndim);
        data.part_stat_list[stat_tag::ener_thm][iker][0] += ener_loss / data.ener_start;

        // compute bin indices
        idx_ener = findIdx(event->ener, bin_list[bin_tag::ener]);

        // compute transport
        dt = event->time - data.time_prev;
        dsplus = event->splus - data.splus_prev;
        dsminus = event->sminus - data.sminus_prev;
        ds = dsplus + dsminus;
        sign = dsplus > dsminus ? 1.0 : -1.0;

        // update B-field scattering
        while ( true ) {
          if ( ds < data.s_scat ) {
            data.pos = data.pos + sign * ds * data.Bhat;
            data.s_scat = data.s_scat - ds;
            break;
          } else {
            data.pos = data.pos + sign * data.s_scat * data.Bhat;
            ds = ds - data.s_scat;
            data.s_scat = -log(1.0 - xi()) * data.lam_scat;
            data.Bhat = calcRandVec(data.mach_A, data.super); // resample B-field direction
          }
        }

        // update time and distance
        data.ener_prev = event->ener;
        data.time_prev = event->time;
        data.splus_prev = event->splus;
        data.sminus_prev = event->sminus;

        // compute energy histograms
        if ( idx_ener > 0 && idx_ener < bin_list[bin_tag::ener].size() ) {
          data.part_stat_list[stat_tag::time_ener][iker][idx_ener - 1] += dt;
        }
        
        // compute interaction histograms
        flag = event->interaction;
        switch ( flag ) {
          case flags::scat: // scattering
          break;
          case flags::brem: // Bremsstrahlung
          data.part_stat_list[stat_tag::ener_loss_mech][iker][flags::brem-1] += event->ener_loss / data.ener_start;
          break;
          case flags::exc: // excitation
          data.part_stat_list[stat_tag::ener_loss_mech][iker][flags::exc-1] += event->ener_loss / data.ener_start;
          break;
          case flags::ion: // ionization
          data.part_stat_list[stat_tag::ener_loss_mech][iker][flags::ion-1] += event->ener_loss / data.ener_start;
          data.part_stat_list[stat_tag::num_ion_elem][iker][event->Zelem - 1] += 1.0;
          // compute secondary energy histograms
          idx_ener_sec = findIdx(event->ener_sec, bin_list[bin_tag::ener_sec]);
          if (idx_ener_sec > 0 && idx_ener_sec < bin_list[bin_tag::ener_sec].size()) {
            data.part_stat_list[stat_tag::ener_thm][iker][0] -= event->ener_sec / data.ener_start; // don't include secondary electron energy in thermalization efficiency
            data.part_stat_list[stat_tag::ener_sec][iker][idx_ener_sec - 1] += event->ener_sec / data.ener_start;
          }
          break;
          case flags::moller: // Moller
          if ( !std::isnan(event->ener_loss) ) { // for some reason, this is sometimes NaN
            data.part_stat_list[stat_tag::ener_loss_mech][iker][flags::moller-1] += event->ener_loss / data.ener_start;
          }
          break;
        }

        // add continuous energy losses
        data.part_stat_list[stat_tag::ener_loss_mech][iker][flags::moller-1] += event->ener_loss_moller / data.ener_start;

        // end conditions
        if ( time > data.dt ) {
          data.outoftime = true;
          flag = flags::outoftime;
        } else if ( event->ener < data.ener_min ) {
          data.thermalized = true;
          flag = flags::thermalized;
        }

        // write data
        if ( do_hist ) {
          data.oss << std::setprecision(15) << time << ",";
          data.oss << std::setprecision(15) << data.pos.x << "," << std::setprecision(15) << data.pos.y << "," << std::setprecision(15) << data.pos.z << ",";
          data.oss << std::setprecision(15) << event->cos_alpha << ", ";
          data.oss << std::setprecision(15) << event->ener << ", ";
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
 * @param M1_stat_list_flat   The flattened list of M1 statistics.
 * @param M2_stat_list_flat   The flattened list of M2 statistics.
 * @param M3_stat_list_flat   The flattened list of M3 statistics.
 * @param M4_stat_list_flat   The flattened list of M4 statistics.
 * @param size_t              The size of the flattened list of statistics.
 */
void getFlatData(
  const vector3d<Data>& data_grid,
  const std::vector<Stat> &stat_list, 
  const size_t nker,
  const size_t nmom,
  std::vector<double> &M1_stat_list_flat, 
  std::vector<double> &M2_stat_list_flat, 
  std::vector<double> &M3_stat_list_flat, 
  std::vector<double> &M4_stat_list_flat, 
  size_t &size_flat
) {
  size_flat = 0;
  size_t size_grid = data_grid.size() * data_grid[0].size() * data_grid[0][0].size();
  for ( size_t ii = 0; ii < stat_list.size(); ii++ ) {
    size_flat += size_grid * nker * stat_list[ii].size;
  }
  M1_stat_list_flat.reserve(size_flat);
  if ( nmom >= 2 ) M2_stat_list_flat.reserve(size_flat);
  if ( nmom >= 3 ) M2_stat_list_flat.reserve(size_flat);
  if ( nmom >= 4 ) M2_stat_list_flat.reserve(size_flat);

  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      for ( size_t k = 0; k < data_grid[i][j].size(); k++ ) {
        const Data &data = data_grid[i][j][k];
        for ( size_t ii = 0; ii < stat_list.size(); ii++ ) {
          const Stat &stat = stat_list[ii];
          for ( size_t jj = 0; jj < nker; jj++ ) {
            for ( size_t kk = 0; kk < stat.size; kk++ ) {
              M1_stat_list_flat.push_back(data.M1_stat_list[ii][jj][kk]);
              if ( nmom >= 2 ) M2_stat_list_flat.push_back(data.M2_stat_list[ii][jj][kk]);
              if ( nmom >= 3 ) M3_stat_list_flat.push_back(data.M3_stat_list[ii][jj][kk]);
              if ( nmom >= 4 ) M4_stat_list_flat.push_back(data.M4_stat_list[ii][jj][kk]);
            }
          }
        }
      }
    }
  }
}
