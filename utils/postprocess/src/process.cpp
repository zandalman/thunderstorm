// includes
#include <fstream>
#include <sstream>
#include <iostream>
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

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A constructor to initialize the Data structure.
Data::Data(double ener_, double escape_, double mach_A_, double L_, int geo_, const std::vector<Stat> &stat_list)
  : ener(ener_)
  , escape(escape_)
  , mach_A(mach_A_)
  , L(L_)
  , geo(geo_)
  , time_start(0.0)
  , splus_start(0.0)
  , sminus_start(0.0)
  , time(0.)
  , splus(0.)
  , sminus(0.)
  , rpar(0.)
  , varpar(0.)
  , varperp(0.)
  , fesc_min(1.0)
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
 }

/// @brief Reset the particle data.
void Data::reset() {
  time_start = 0.;
  splus_start = 0.;
  sminus_start = 0.;
  time = 0.;
  splus = 0.;
  sminus = 0.;
  rpar = 0.;
  varpar = 0.;
  varperp = 0.;
  fesc_min = 1.0;
  for ( size_t i = 0; i < part_stat_list.size(); i++ ) {
    std::fill(part_stat_list[i].begin(), part_stat_list[i].end(), 0.0);
  }
}

void Data::calcStat(int n_int, const std::vector<Stat> &stat_list) {
  Stat stat;
  double delta, delta_np1, delta_np1_sq, term1;
  double n = static_cast<double>(n_int);
  double np1 = n + 1.0;
  for ( size_t i = 0; i < stat_list.size(); i++ ) {
    stat = stat_list[i];
    for ( size_t j = 0; j < stat.size; j++ ) {
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

/**
 * @brief Process a binary data file output by thunderstorm.
 * 
 * @param datafile_name       The name of the binary data file.
 * @param num_event_per_chunk The number of events per chunk.
 * @param count               The particle count.
 * @param data_grid           The grid of data.  
 */
void processFile(const std::string &datafile_name, size_t num_event_per_chunk, const vector2d<double> &bin_list, const std::vector<Stat>& stat_list, int &count, vector2d<Data>& data_grid) {

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
  while ( datafile.read(buffer.data(), chunk_size) ) {
    for ( size_t i = 0; i < num_event_per_chunk; i++ ) {
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if ( event->id != current_id ) {
        count++;
        for ( size_t j = 0; j < data_grid.size(); j++ ) {
          for ( size_t k = 0; k < data_grid[j].size(); k++ ) {
            data_grid[j][k].calcStat(count, stat_list);
            data_grid[j][k].reset();
          }
        }
        current_id = event->id;
      } else {
        processEvent(event, bin_list, stat_list, data_grid);
      }
    }
  }

  // Handle the case where the last chunk might not be full
  if ( datafile.eof() ) {
    size_t bytes_read = datafile.gcount();
    if ( bytes_read > 0 ) {
      size_t num_event_last_chunk = bytes_read / event_size;
      for ( size_t i = 0; i < num_event_last_chunk; i++ ) {
        event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
        if ( event->id != current_id ) {
          count++;
          for ( size_t j = 0; j < data_grid.size(); j++ ) {
            for ( size_t k = 0; k < data_grid[j].size(); k++ ) {
              data_grid[j][k].calcStat(count, stat_list);
              data_grid[j][k].reset();
            }
          }
          current_id = event->id;
        } else {
          processEvent(event, bin_list, stat_list, data_grid);
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
 * @brief Process an event and add the post-processed data to the particle data.
 * 
 * @param event     The event struct.
 * @param data_grid The grid of data.
 */
void processEvent(const Event* event, const vector2d<double> &bin_list, const std::vector<Stat>& stat_list, vector2d<Data>& data_grid) {

  if ( event == nullptr ) {
    std::cerr << "Error: Null event pointer." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  double time_rel, splus_rel, sminus_rel, s_rel;
  double dt, dsplus, dsminus;
  double fesc, ener_loss;
  size_t idx_ener_sec, idx_time;
  
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
    
      Data &data = data_grid[i][j];
      if ( data.fesc_min == 0.0 ) continue;
      
      // if energy is below starting energy, update start time and coordinates
      if ( event->ener > data.ener ) {
        data.time_start = event->time;
        data.splus_start = event->splus;
        data.sminus_start = event->sminus;
        data.time = event->time;
        data.splus = event->splus;
        data.sminus = event->sminus;
        continue;
      }

      // compute useful quantities
      time_rel = event->time - data.time_start;
      splus_rel = event->splus - data.splus_start;
      sminus_rel = event->sminus - data.sminus_start;
      s_rel = splus_rel + sminus_rel;

      ener_loss = event->ener_loss_moller + event->ener_loss_sync + event->ener_loss_cher;
      if ( !std::isnan(event->ener_loss) ) {
        ener_loss += event->ener_loss;
      }

      // compute transport
      dt = event->time - data.time;
      dsplus = event->splus - data.splus;
      dsminus = event->sminus - data.sminus;
      calcTransport(data.mach_A, data.L, s_rel, dsplus, dsminus, data.rpar, data.varpar, data.varperp);
      calcFesc(data.geo, data.escape, data.rpar, data.varpar, data.varperp, fesc);
      if ( fesc < data.fesc_min ) {
        data.fesc_min = fesc;
      } else {
        fesc = data.fesc_min;
      }

      // update time and distance
      data.time = event->time;
      data.splus = event->splus;
      data.sminus = event->sminus;

      // compute time histograms
      idx_time = findIdx(time_rel, bin_list[bin_tag::time]);
      if ( idx_time > 0 && idx_time < bin_list[bin_tag::time].size() ) {
        data.part_stat_list[stat_tag::ener_loss_time][idx_time - 1] += fesc * ener_loss;
      }
      
      // compute interaction histograms
      switch ( event->interaction ) {
        case flags::brem:
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::brem] += fesc * event->ener_loss;
        break;
        case flags::exc:
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::exc] += fesc * event->ener_loss;
        break;
        case flags::ion:
        data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::ion] += fesc * event->ener_loss;
        data.part_stat_list[stat_tag::num_ion_elem][event->Zelem - 1] += fesc;
        // compute secondary energy histograms
        idx_ener_sec = findIdx(event->ener_sec, bin_list[bin_tag::ener_sec]);
        if (idx_ener_sec > 0 && idx_ener_sec < bin_list[bin_tag::ener_sec].size()) {
          data.part_stat_list[stat_tag::num_sec_ener][idx_ener_sec - 1] += fesc;
        }
        break;
        case flags::moller:
        if ( !std::isnan(event->ener_loss) ) {
          data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::moller] += fesc * event->ener_loss;
        }
        break;
      }
      data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::moller] += fesc * event->ener_loss_moller;
      data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::sync] += fesc * event->ener_loss_sync;
      data.part_stat_list[stat_tag::ener_loss_mech][mech_tag::cher] += fesc * event->ener_loss_cher;
    }
  }
}

/**
 * @brief Flatten the data so it can be communicated via MPI.
 * 
 * @param data_grid           The grid of data.
 * @param ener_loss_mech_flat The energy loss for each mechanism [eV].
 * @param num_ion_elem        The number of ionizations per element.
 * @param num_sec_ener        The number of secondary particles per energy bin.
 */
void getFlatData(const vector2d<Data>& data_grid, const std::vector<Stat> &stat_list, std::vector<double> &mean_stat_list_flat, std::vector<double> &M2_stat_list_flat, std::vector<double> &M3_stat_list_flat, std::vector<double> &M4_stat_list_flat, size_t &size_flat) {
  Data data;
  Stat stat;
  size_flat = 0;
  for ( size_t i = 0; i < data_grid.size(); i++ ) {
    for ( size_t j = 0; j < data_grid[i].size(); j++ ) {
      data = data_grid[i][j];
      for ( size_t k = 0; k < stat_list.size(); k++ ) {
        stat = stat_list[k];
        for ( size_t l = 0; l < stat.size; l++ ) {
          mean_stat_list_flat.push_back(data.mean_stat_list[k][l]);
          M2_stat_list_flat.push_back(data.M2_stat_list[k][l]);
          M3_stat_list_flat.push_back(data.M2_stat_list[k][l]);
          M4_stat_list_flat.push_back(data.M2_stat_list[k][l]);
          size_flat += 1;
        }
      }
    }
  }
}
