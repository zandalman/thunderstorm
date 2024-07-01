// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>

// headers
#include "const.h"
#include "functions.h"
#include "io.h"

/// @brief A constructor to initialize the Event structure.
Event::Event(int int_data[5], double double_data[11]) 
  : id(int_data[0])                   // The particle ID.
  , nstep(int_data[1])                // The step number.
  , Zelem(int_data[2])                // The proton number of the element.
  , interaction(int_data[3])          // The interaction flag.
  , ion(int_data[4])                  // The ion index.
  , time(double_data[0])              // The event time [s].
  , x(double_data[1])                 // The event x-coordinate [cm].
  , y(double_data[2])                 // The event y-coordinate [cm].
  , z(double_data[3])                 // The event z-coordinate [cm].
  , ener(double_data[4])              // The particle kinetic energy [eV].
  , cos_th(double_data[5])            // The cosine of the scattering angle.
  , ener_loss(double_data[6])         // The energy lost [eV].
  , ener_sec(double_data[7])          // The energy of the secondary [eV].
  , ener_loss_sync(double_data[8])    // The energy lost due to synchrotron [eV].
  , ener_loss_cher(double_data[9])    // The energy lost due to Cherenkov radiation [eV].
  , ener_loss_moller(double_data[10]) // The energy lost due to small-angle Moller scattering [eV].
 {}

PartData::PartData(size_t num_time, size_t num_dis)
  : id(-1)
 {
  for ( size_t i = 0; i < 6; i++ ) ener_loss_mech[i] = 0.;
  for ( size_t i = 0; i < 118; i++ ) num_ion[i] = 0;
  for ( size_t i = 0; i < (num_time - 1); i++ ) num_ev_time.push_back(0);
  for ( size_t i = 0; i < (num_time - 1); i++ ) ener_time.push_back(0.);
  for ( size_t i = 0; i < (num_time - 1); i++ ) ener_loss_time.push_back(0.);
  for ( size_t i = 0; i < (num_dis - 1); i++ ) ener_loss_dis.push_back(0.);
 }

void PartData::reset(int id_, double ener_, size_t num_time, size_t num_dis) {
  id = id_;
  num_ev = 0;
  ener = ener_;
  x_start = 0.;
  y_start = 0.;
  z_start = 0.;
  for ( size_t i = 0; i < 6; i++ ) ener_loss_mech[i] = 0.;
  for ( size_t i = 0; i < 118; i++ ) num_ion[i] = 0;
  for ( size_t i = 0; i < (num_time - 1); i++ ) num_ev_time[i] = 0;
  for ( size_t i = 0; i < (num_time - 1); i++ ) ener_time[i] = 0.;
  for ( size_t i = 0; i < (num_time - 1); i++ ) ener_loss_time[i] = 0.;
  for ( size_t i = 0; i < (num_dis - 1); i++ ) ener_loss_dis[i] = 0.;
}

void processEvent(const Event* event, std::vector<PartData>& part_data_list, const std::vector<double> &ener_list, const std::vector<double> &time_list, const std::vector<double> &dis_list) {

  if ( event == nullptr ) {
    std::cerr << "Error: Null event pointer." << std::endl;
    return;
  }

  for ( size_t i = 0; i < ener_list.size(); i++ ) {
    
    PartData &part_data = part_data_list[i];
    
    if ( event->ener > ener_list[i] ) {
      part_data.x_start = event->x;
      part_data.y_start = event->y;
      part_data.z_start = event->z;
      continue;
    }

    // increment the event count
    part_data.num_ev += 1;
    
    // compute time index
    size_t idx_time = findIdx(event->time, time_list);
    
    // compute distance index
    double x_rel = event->x - part_data.x_start;
    double y_rel = event->y - part_data.y_start;
    double z_rel = event->z - part_data.z_start;
    double dis = sqrt(x_rel*x_rel + y_rel*y_rel + z_rel*z_rel);
    size_t idx_dis = findIdx(dis, dis_list);

    // compute energy loss
    double ener_loss = event->ener_loss + event->ener_loss_moller + event->ener_loss_sync + event->ener_loss_cher;

    // compute time histograms
    if ( idx_time > 0 && idx_time < time_list.size() ) {
      part_data.num_ev_time[idx_time] += 1;
      part_data.ener_time[idx_time] += event->ener;
      part_data.ener_loss_time[idx_time] += ener_loss;
    }

    // compute distance histograms
    if ( idx_dis > 0 && idx_dis < dis_list.size() ) {
      part_data.ener_loss_dis[idx_dis] += ener_loss;
    }
    
    // compute interaction histograms
    switch ( event->interaction ) {
      case flags::brem:
      part_data.ener_loss_mech[0] += event->ener_loss;
      break;
      case flags::exc:
      part_data.ener_loss_mech[1] += event->ener_loss;
      break;
      case flags::ion:
      part_data.ener_loss_mech[2] += event->ener_loss;
      part_data.num_ion[event->Zelem - 1] += 1;
      break;
      case flags::moller:
      part_data.ener_loss_mech[3] += event->ener_loss;
      break;
    }
    part_data.ener_loss_mech[3] += event->ener_loss_moller;
    part_data.ener_loss_mech[4] += event->ener_loss_sync;
    part_data.ener_loss_mech[5] += event->ener_loss_cher;

  }
}

void writePartData(std::vector<PartData>& part_data_list, const std::vector<double> &ener_list, const std::vector<double> &time_list, const std::vector<double> &dis_list) {
  for ( size_t i = 0; i < ener_list.size(); i++ ) {
    PartData &part_data = part_data_list[i];
    for ( size_t j = 0; j < (time_list.size() - 1); j++ ) {
      if ( part_data.num_ev_time[j] > 0 ) {
        part_data.ener_time[j] /= part_data.num_ev_time[j];
      }
    }
  }
}

int processFile(std::string filename, size_t num_event_per_chunk, const std::vector<double> &ener_list, const std::vector<double> &time_list, const std::vector<double> &dis_list) {

  // Open file
  std::ifstream file(filename, std::ios::binary);

  if ( !file.is_open() ) {
    std::cerr << "Failed to open file " << filename << std::endl;
    return 1;
  }

  // Compute chunk size
  const size_t event_size = sizeof(Event);
  const size_t chunk_size = event_size * num_event_per_chunk;

  // Allocate memory to hold the chunk content
  std::vector<char> buffer(chunk_size);
  Event* event;

  size_t npart = 0;
  int current_id = -1;
  std::vector<PartData> part_data_list(ener_list.size(), PartData(time_list.size(), dis_list.size()));
  
  while ( file.read(buffer.data(), chunk_size) ) {
    for ( size_t i = 0; i < num_event_per_chunk; i++ ) {
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if ( event->id != current_id ) {
        if ( current_id != -1 ) writePartData(part_data_list, ener_list, time_list, dis_list);
        for ( size_t i = 0; i < ener_list.size(); i++ ) part_data_list[i].reset(event->id, ener_list[i], time_list.size(), dis_list.size());
        current_id = event->id;
        npart++;
      } else {
        processEvent(event, part_data_list, ener_list, time_list, dis_list);
      }
    }
  }

  if ( file.eof() ) {
    size_t bytes_read = file.gcount();
    // Handle the case where the last chunk might not be full
    if ( bytes_read > 0 ) {
      size_t num_event_last_chunk = bytes_read / event_size;
      for ( size_t i = 0; i < num_event_last_chunk; i++ ) {
        event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
        if ( event->id != current_id ) {
          if ( current_id != -1 ) writePartData(part_data_list, ener_list, time_list, dis_list);
          for ( size_t i = 0; i < ener_list.size(); i++ ) part_data_list[i].reset(event->id, ener_list[i], time_list.size(), dis_list.size());
          current_id = event->id;
          npart++;
        } else {
          processEvent(event, part_data_list, ener_list, time_list, dis_list);
        }
      }
    }
  } else {
    std::cerr << "Error reading file " << filename << std::endl;
    return 1;
  }

  file.close();
  return 0;
}
