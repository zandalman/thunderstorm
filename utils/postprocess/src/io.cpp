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
#include "io.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

/// @brief A constructor to initialize the Event structure.
Event::Event(int int_data[5], double double_data[12]) 
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
  , cos_alpha(double_data[6])         // The particle pitch angle cosine.
  , ener_loss(double_data[7])         // The energy lost [eV].
  , ener_sec(double_data[8])          // The energy of the secondary [eV].
  , ener_loss_sync(double_data[9])    // The energy lost due to synchrotron [eV].
  , ener_loss_cher(double_data[10])    // The energy lost due to Cherenkov radiation [eV].
  , ener_loss_moller(double_data[11]) // The energy lost due to small-angle Moller scattering [eV].
 {}

/// @brief A constructor to initialize the PartData structure.
PartData::PartData(size_t num_ener_sec, size_t num_time, size_t num_dis_par, size_t num_dis_perp, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat)
  : id(-1)
 {
  for ( size_t i = 0; i < 6; i++ ) {
    ener_loss_mech.push_back(0.);
  }
  for ( size_t i = 0; i < 118; i++ ) {
    num_ion_elem.push_back(0);
  }
  if ( do_stat_cat[0] ) {
    for ( size_t i = 0; i < (num_ener_sec - 1); i++ ) {
      if ( do_stat_list[2] ) num_sec_ener.push_back(0);
    }
  }
  if ( do_stat_cat[1] ) {
    for ( size_t i = 0; i < (num_time - 1); i++ ) {
      num_ev_time.push_back(0);
      if ( do_stat_list[3] ) dis_par_time.push_back(0.);
      if ( do_stat_list[4] ) dis_perp_time.push_back(0.);
      if ( do_stat_list[5] ) ener_time.push_back(0.);
      if ( do_stat_list[6] ) ener_loss_time.push_back(0.);
    }
  }
  if ( do_stat_cat[2] && do_stat_cat[3] ) {
    std::vector<double> ener_loss_dis_perp;
    for ( size_t i = 0; i < (num_dis_par - 1); i++ ) {
      ener_loss_dis_perp.clear();
      for ( size_t j = 0; j < (num_dis_perp - 1); j++ ) {
        if ( do_stat_list[7] ) ener_loss_dis_perp.push_back(0.);
      }
      if ( do_stat_list[7] ) ener_loss_dis2d.push_back(ener_loss_dis_perp);
    }
  }
 }

/**
 * @brief Reset the particle data.
 * 
 * @param id_          The new particle ID.
 * @param ener_        The new particle energy.
 * @param num_ener_sec The number of secondary energy bins.
 * @param num_time     The number of time bins.
 * @param num_dis      The number of distance bins.
 */
void PartData::reset(int id_, double ener_, size_t num_ener_sec, size_t num_time, size_t num_dis_par, size_t num_dis_perp, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat) {
  id = id_;
  num_ev = 0;
  ener = ener_;
  t_start = 0.;
  x_start = 0.;
  y_start = 0.;
  z_start = 0.;
  for ( size_t i = 0; i < 6; i++ ) {
    ener_loss_mech[i] = 0.;
  }
  for ( size_t i = 0; i < 118; i++ ) {
    num_ion_elem[i] = 0;
  }
  if ( do_stat_cat[0] ) {
    for ( size_t i = 0; i < (num_ener_sec - 1); i++ ) {
      if ( do_stat_list[2] ) num_sec_ener[i] = 0;
    }
  }
  if ( do_stat_cat[1] ) {
    for ( size_t i = 0; i < (num_time - 1); i++ ) {
      num_ev_time[i] = 0;
      if ( do_stat_list[3] ) dis_par_time[i] = 0.;
      if ( do_stat_list[4] ) dis_perp_time[i] = 0.;
      if ( do_stat_list[5] ) ener_time[i] = 0.;
      if ( do_stat_list[6] ) ener_loss_time[i] = 0.;
    }
  }
  if ( do_stat_cat[2] && do_stat_cat[3] ) {
    for ( size_t i = 0; i < (num_dis_par - 1); i++ ) {
      for ( size_t j = 0; j < (num_dis_perp - 1); j++ ) {
        if ( do_stat_list[7] ) ener_loss_dis2d[i][j] = 0.;
      }
    }
  }
}

/**
 * @brief Clear a file.
 * 
 * @param file_name The name of the file to clear.
*/
void clearFile(const std::string& file_name) {
  std::ofstream file(file_name);
  file.close();
}

/**
 * @brief Write a vector to a string stream.
 * 
 * @param oss The string stream.
 * @param vec The vector to write.
 */
template <typename T>
void writeVector(std::ostringstream &oss, const std::vector<T>& vec) {
  for (size_t i = 0; i < vec.size(); i++) {
    oss << vec[i];
    if (i != vec.size() - 1) oss << ",";
  }
  oss << std::endl;
}

/**
 * @brief Write a 2d vector to a file.
 * 
 * @param oss The string stream.
 * @param vec The vector to write.
 */
template <typename T>
void writeVector2d(std::ostringstream &oss, const vector2d<T>& vec) {
  for (size_t i = 0; i < vec.size(); ++i) {
    for (size_t j = 0; j < vec[i].size(); j++) {
      oss << vec[i][j];
      if (j != vec[i].size() - 1) oss << ",";
    }
    if (i != vec.size() - 1) oss << ";";
  }
  oss << std::endl;
}

/**
 * @brief Write information about the post-processing to a text file.
 * 
 * @param infofile_name The name of the information file.
 * @param config        The struct containing the configuration file information.
 * @param ener_list     The list of energy bins.
 * @param ener_sec_list The list of secondary energy bins.
 * @param time_list     The list of time bins.
 * @param dis_list      The list of distance bins.
 */
void writeInfo(std::string &infofile_name, Config &config, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list, size_t num_stat, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat) {

  std::ostringstream oss;

  oss << "Thunderstorm post-processing" << std::endl;
  oss << std::endl;

  oss << "Post-processing parameters" << std::endl;
  oss << "File number:        " << config["IO"]["num_file"] << std::endl;
  oss << "Events per chunk:   " << config["IO"]["num_event_per_chunk"] << std::endl;
  oss << "Number of energies: " << config["Energy"]["num"] << std::endl;
  oss << "Number of lines:    " << (num_stat + 4) << std::endl;
  oss << std::endl;

  oss << "Variable lists" << std::endl;
  oss << "Energy [eV]" << std::endl;  
  writeVector(oss, ener_list);
  if ( do_stat_cat[0] ) {
    oss << "Secondary energy [eV]" << std::endl;  
    writeVector(oss, ener_sec_list);
  }
  if ( do_stat_cat[1] ) {
    oss << "Time [s]" << std::endl;
    writeVector(oss, time_list);
  }
  if ( do_stat_cat[2] ) {
    oss << "Distance (parallel) [cm]" << std::endl;
    writeVector(oss, dis_par_list);
  }
  if ( do_stat_cat[3] ) {
    oss << "Distance (perppendicular) [cm]" << std::endl;
    writeVector(oss, dis_perp_list);
  }
  oss << std::endl;

  size_t num_line = 1;
  oss << "Post-processed data format" << std::endl;
  oss << "1.  energy" << std::endl;
  oss << "2.  number of events" << std::endl;
  oss << "3.  start time [s] and start coordinates [cm]" << std::endl;
  oss << "4.  death time [s] and death coordinates [cm]" << std::endl;
  num_line = 5;

  if ( do_stat_list[0] ) oss << num_line << ".  " << "energy loss [eV] for each mechanism" << std::endl, num_line++;
  if ( do_stat_list[1] ) oss << num_line << ".  " << "number of ionizations per element" << std::endl, num_line++;
  if ( do_stat_list[2] ) oss << num_line << ".  " << "number of secondary electrons per secondary energy bin" << std::endl, num_line++;
  if ( do_stat_list[3] ) oss << num_line << ".  " << "average parallel distance [cm] for each time bin" << std::endl, num_line++;
  if ( do_stat_list[4] ) oss << num_line << ".  " << "average perpendicular distance [cm] for each time bin" << std::endl, num_line++;
  if ( do_stat_list[5] ) oss << num_line << ".  " << "average energy [eV] for each time bin" << std::endl, num_line++;
  if ( do_stat_list[6] ) oss << num_line << ".  " << "energy loss [eV] per time bin" << std::endl, num_line++;
  if ( do_stat_list[7] ) oss << num_line << ".  " << "energy loss [eV] per parallel and perpendicular distance bin" << std::endl, num_line++;

  std::ofstream infofile(infofile_name);
  if ( !infofile ) {
    std::cerr << "Failed to open file " << infofile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  infofile << oss.str();
  infofile.close();
  if ( !infofile.good() ) {
    std::cerr << "Error writing to file " << infofile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

/**
 * @brief Process an event and add the post-processed data to the particle data.
 * 
 * @param event          The event struct.
 * @param part_data_list The list of particle data structs for each energy.
 * @param ener_list      The list of energy bins.
 * @param ener_sec_list  The list of secondary energy bins.
 * @param time_list      The list of time bins.
 * @param dis_list       The list of distance bins.
 */
void processEvent(const Event* event, std::vector<PartData>& part_data_list, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat) {

  if ( event == nullptr ) {
    std::cerr << "Error: Null event pointer." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  double t_rel, x_rel, y_rel, z_rel, dis_perp, ener_loss;
  size_t idx_ener_sec, idx_time, idx_dis_par, idx_dis_perp;
  
  for ( size_t i = 0; i < ener_list.size(); i++ ) {
    
    PartData &part_data = part_data_list[i];
    
    // if energy is below starting energy, update start time and coordinates
    if ( event->ener > ener_list[i] ) {
      part_data.t_start = event->time;
      part_data.x_start = event->x;
      part_data.y_start = event->y;
      part_data.z_start = event->z;
      continue;
    }
    
    // compute useful quantities
    part_data.num_ev += 1;
    t_rel = event->time - part_data.t_start;
    x_rel = event->x - part_data.x_start;
    y_rel = event->x - part_data.y_start;
    z_rel = event->x - part_data.z_start;
    dis_perp = sqrt(x_rel*x_rel + y_rel*y_rel);
    ener_loss = event->ener_loss + event->ener_loss_moller + event->ener_loss_sync + event->ener_loss_cher;

    // compute indices
    if ( do_stat_cat[1] ) idx_time = findIdx(t_rel, time_list);
    if ( do_stat_cat[2] ) idx_dis_par = findIdx(z_rel, dis_par_list);
    if ( do_stat_cat[3] ) idx_dis_perp = findIdx(dis_perp, dis_perp_list);

    // compute time histograms
    if ( do_stat_cat[1] && (idx_time > 0) && (idx_time < time_list.size()) ) {
      part_data.num_ev_time[idx_time - 1] += 1;
      if ( do_stat_list[3] ) part_data.dis_par_time[idx_time - 1] += z_rel;
      if ( do_stat_list[4] ) part_data.dis_perp_time[idx_time - 1] += dis_perp;
      if ( do_stat_list[5] ) part_data.ener_time[idx_time - 1] += event->ener;
      if ( do_stat_list[6] ) part_data.ener_loss_time[idx_time - 1] += ener_loss;
    }

    // compute 2d distance histograms
    if ( do_stat_cat[2] && (idx_dis_par > 0) && (idx_dis_par < dis_par_list.size()) && (idx_dis_perp > 0) && (idx_dis_par < dis_perp_list.size())) {
      if ( do_stat_list[7] ) part_data.ener_loss_dis2d[idx_dis_par - 1][idx_dis_perp - 1] += ener_loss;
    }
    
    // compute interaction histograms
    switch ( event->interaction ) {
      case flags::death:
      part_data.t_end = t_rel;
      part_data.x_end = x_rel;
      part_data.y_end = y_rel;
      part_data.z_end = z_rel;
      break;
      case flags::brem:
      if ( do_stat_list[0] ) part_data.ener_loss_mech[0] += event->ener_loss;
      break;
      case flags::exc:
      if ( do_stat_list[0] ) part_data.ener_loss_mech[1] += event->ener_loss;
      break;
      case flags::ion:
      if ( do_stat_list[0] ) part_data.ener_loss_mech[2] += event->ener_loss;
      if ( do_stat_list[1] ) part_data.num_ion_elem[event->Zelem - 1] += 1;
      
      // compute secondary energy histograms
      if ( do_stat_cat[0] ) {
        idx_ener_sec = findIdx(event->ener_sec, ener_sec_list);
        if (idx_ener_sec > 0 && idx_ener_sec < ener_sec_list.size()) {
          if ( do_stat_list[2] ) part_data.num_sec_ener[idx_ener_sec - 1] += 1;
        }
      }

      break;
      case flags::moller:
      if ( do_stat_list[0] ) part_data.ener_loss_mech[3] += event->ener_loss;
      break;
    }
    if ( do_stat_list[0] ) {
      part_data.ener_loss_mech[3] += event->ener_loss_moller;
      part_data.ener_loss_mech[4] += event->ener_loss_sync;
      part_data.ener_loss_mech[5] += event->ener_loss_cher;
    }
  }
}

/**
 * @brief Post-process particle data.
 * 
 * @param part_data_list The list of particle data structs for each energy.
 */
void postProcPartData(std::vector<PartData>& part_data_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat) {
  if ( do_stat_cat[1] ) {
    for ( size_t i = 0; i < part_data_list.size(); i++ ) {
      PartData &part_data = part_data_list[i];
      for ( size_t j = 0; j < (part_data.num_ev_time.size() - 1); j++ ) {
        if ( part_data.num_ev_time[j] > 0 ) {
          if ( do_stat_list[3] ) part_data.dis_par_time[j] /= part_data.num_ev_time[j];
          if ( do_stat_list[4] ) part_data.dis_perp_time[j] /= part_data.num_ev_time[j];
          if ( do_stat_list[5] ) part_data.ener_time[j] /= part_data.num_ev_time[j];
        } else {
          if ( do_stat_list[3] ) part_data.dis_par_time[j] = -1.;
          if ( do_stat_list[4] ) part_data.dis_perp_time[j] = -1.;
          if ( do_stat_list[5] ) part_data.ener_time[j] = -1.;
        }
      }
    }
  }
}

/**
 * @brief Write particle data to an output file.
 * 
 * @param outfile_name   The name of the output file.
 * @param part_data_list The list of particle data structs for each energy.
 */
void writePartData(std::string &outfile_name, std::vector<PartData>& part_data_list, const std::vector<bool> &do_stat_list) {

  std::ostringstream oss;

  for ( size_t i = 0; i < part_data_list.size(); i++ ) {
    PartData &part_data = part_data_list[i];

    oss << part_data.ener << std::endl;
    oss << part_data.num_ev << std::endl;
    oss << part_data.t_start << ", " << part_data.x_start << "," << part_data.y_start << "," << part_data.z_start << std::endl;
    oss << part_data.t_end << ", " << part_data.x_end << "," << part_data.y_end << "," << part_data.z_end << std::endl;
    if ( do_stat_list[0] ) writeVector(oss, part_data.ener_loss_mech);
    if ( do_stat_list[1] ) writeVector(oss, part_data.num_ion_elem);
    if ( do_stat_list[2] ) writeVector(oss, part_data.num_sec_ener);
    if ( do_stat_list[3] ) writeVector(oss, part_data.dis_par_time);
    if ( do_stat_list[4] ) writeVector(oss, part_data.dis_perp_time);
    if ( do_stat_list[5] ) writeVector(oss, part_data.ener_time);
    if ( do_stat_list[6] ) writeVector(oss, part_data.ener_loss_time);
    if ( do_stat_list[7] ) writeVector2d(oss, part_data.ener_loss_dis2d);
  }

  std::ofstream outfile(outfile_name, std::ios::app);

  if ( !outfile ) {
    std::cerr << "Failed to open file for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  outfile << oss.str();
  outfile.close();
  if ( !outfile.good() ) {
    std::cerr << "Error writing to file." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

/**
 * @brief Process a binary data file output by thunderstorm.
 * 
 * @param datafile_name       The name of the binary data file.
 * @param outfile_name        The name of the output file.
 * @param num_event_per_chunk The number of events per chunk.
 * @param ener_list           The list of energy bins.
 * @param ener_sec_list       The list of secondary energy bins.
 * @param time_list           The list of time bins.
 * @param dis_list            The list of distance bins.    
 */
int processFile(std::string& datafile_name, std::string& outfile_name, size_t num_event_per_chunk, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list, const std::vector<bool> &do_stat_list, const std::vector<bool> &do_stat_cat) {

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

  int npart = 0;
  int current_id = -1;
  std::vector<PartData> part_data_list(ener_list.size(), PartData(ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size(), do_stat_list, do_stat_cat));
  
  while ( datafile.read(buffer.data(), chunk_size) ) {
    for ( size_t i = 0; i < num_event_per_chunk; i++ ) {
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if ( event->id != current_id ) {
        if ( current_id != -1 ) {
          postProcPartData(part_data_list, do_stat_list, do_stat_cat);
          writePartData(outfile_name, part_data_list, do_stat_list);
        }
        for ( size_t j = 0; j < ener_list.size(); j++ ) {
          part_data_list[j].reset(event->id, ener_list[j], ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size(), do_stat_list, do_stat_cat);
        }
        current_id = event->id;
        npart++;
      } else {
        processEvent(event, part_data_list, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list, do_stat_list, do_stat_cat);
      }
    }
  }

  if ( datafile.eof() ) {
    size_t bytes_read = datafile.gcount();
    // Handle the case where the last chunk might not be full
    if ( bytes_read > 0 ) {
      size_t num_event_last_chunk = bytes_read / event_size;
      for ( size_t i = 0; i < num_event_last_chunk; i++ ) {
        event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
        if ( event->id != current_id ) {
          if ( current_id != -1 ) {
            postProcPartData(part_data_list, do_stat_list, do_stat_cat);
            writePartData(outfile_name, part_data_list, do_stat_list);
          }
          for ( size_t j = 0; j < ener_list.size(); j++ ) {
            part_data_list[j].reset(event->id, ener_list[j], ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size(), do_stat_list, do_stat_cat);
          }
          current_id = event->id;
          npart++;
        } else {
          processEvent(event, part_data_list, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list, do_stat_list, do_stat_cat);
        }
      }
    }
  } else {
    std::cerr << "Error reading file " << datafile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  datafile.close();
  return npart;
}

/**
 * @brief Concatenate all output files into a single file.
 * 
 * @param out_path The path of the output file directory.
 * @param num_file The number of output files.
 */
void concatenateFiles(std::string out_path, int num_file) {

  std::string outfile_name = out_path + "/data.txt";
  std::ofstream outfile(outfile_name, std::ios::out | std::ios::trunc);

  if ( !outfile ) {
    std::cerr << "Failed to open " << outfile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for ( int i = 0; i < num_file; i++ ) {
    
    std::string datafile_name = out_path + "/data.txt." + std::to_string(i);
    std::ifstream datafile(datafile_name, std::ios::in);

    if ( !datafile ) {
    std::cerr << "Failed to open " << datafile_name << " for reading." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    outfile << datafile.rdbuf();

    datafile.close();
    if ( !datafile.good() ) {
      std::cerr << "Error reading file " << datafile_name << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if ( std::remove(datafile_name.c_str()) != 0 ) {
      std::cerr << "Error deleting file " << datafile_name << std::endl;
    }
  }
  
  if ( !outfile.good() ) {
    std::cerr << "Error writing to file " << outfile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  outfile.close();
}
