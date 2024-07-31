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
PartData::PartData(size_t num_ener_sec, size_t num_time, size_t num_dis_par, size_t num_dis_perp)
  : id(-1)
 {
  for ( size_t i = 0; i < 6; i++ ) {
    ener_loss_mech.push_back(0.);
  }
  for ( size_t i = 0; i < 118; i++ ) {
    num_ion_elem.push_back(0);
  }
  for ( size_t i = 0; i < (num_ener_sec - 1); i++ ) {
    num_sec_ener.push_back(0);
  }
  for ( size_t i = 0; i < (num_time - 1); i++ ) {
    num_ev_time.push_back(0);
    dis_par_time.push_back(0.);
    dis_perp_time.push_back(0.);
    ener_time.push_back(0.);
    ener_loss_time.push_back(0.);
  }
  std::vector<double> ener_loss_dis_perp;
  for ( size_t i = 0; i < (num_dis_par - 1); i++ ) {
    ener_loss_dis_perp.clear();
    for ( size_t j = 0; j < (num_dis_perp - 1); j++ ) {
      ener_loss_dis_perp.push_back(0.);
    }
    ener_loss_dis2d.push_back(ener_loss_dis_perp);
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
void PartData::reset(int id_, double ener_, size_t num_ener_sec, size_t num_time, size_t num_dis_par, size_t num_dis_perp) {
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
  for ( size_t i = 0; i < (num_ener_sec - 1); i++ ) {
    num_sec_ener[i] = 0;
  }
  for ( size_t i = 0; i < (num_time - 1); i++ ) {
    num_ev_time[i] = 0;
    dis_par_time[i] = 0.;
    dis_perp_time[i] = 0.;
    ener_time[i] = 0.;
    ener_loss_time[i] = 0.;
  }
  for ( size_t i = 0; i < (num_dis_par - 1); i++ ) {
    for ( size_t j = 0; j < (num_dis_perp - 1); j++ ) {
      ener_loss_dis2d[i][j] = 0.;
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
 * @brief Write a vector to a file.
 * 
 * @param file The open file object.
 * @param vec  The vector to write.
 */
template <typename T>
void writeVector(std::ofstream& file, const std::vector<T>& vec) {

  if ( !file.is_open() ) {
    std::cerr << "File is not open." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (size_t i = 0; i < vec.size(); ++i) {
    file << vec[i];
    if (i != vec.size() - 1) file << ",";
  }
  file << std::endl;
}

/**
 * @brief Write a 2d vector to a file.
 * 
 * @param file The open file object.
 * @param vec  The vector to write.
 */
template <typename T>
void writeVector2d(std::ofstream& file, const vector2d<T>& vec) {

  if ( !file.is_open() ) {
    std::cerr << "File is not open." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for (size_t i = 0; i < vec.size(); ++i) {
    for (size_t j = 0; j < vec[i].size(); j++) {
      file << vec[i][j];
      if (j != vec[i].size() - 1) file << ",";
    }
    if (i != vec.size() - 1) file << ";";
  }
  file << std::endl;
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
void writeInfo(std::string &infofile_name, Config &config, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list) {

  std::ofstream infofile(infofile_name);

  if ( !infofile ) {
    std::cerr << "Failed to open file " << infofile_name << " for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  infofile << "Thunderstorm post-processing" << std::endl;
  infofile << std::endl;

  infofile << "Post-processing parameters" << std::endl;
  infofile << "File number:        " << config["IO"]["num_file"] << std::endl;
  infofile << "Events per chunk:   " << config["IO"]["num_event_per_chunk"] << std::endl;
  infofile << "Number of energies: " << config["Energy"]["num"] << std::endl;
  infofile << "Number of lines:    " << 11 << std::endl;
  infofile << std::endl;

  infofile << "Variable lists" << std::endl;
  infofile << "Energy [eV]" << std::endl;  
  writeVector(infofile, ener_list);
  infofile << "Secondary energy [eV]" << std::endl;  
  writeVector(infofile, ener_sec_list);
  infofile << "Time [s]" << std::endl;
  writeVector(infofile, time_list);
  infofile << "Distance (parallel) [cm]" << std::endl;
  writeVector(infofile, dis_par_list);
  infofile << "Distance (perppendicular) [cm]" << std::endl;
  writeVector(infofile, dis_perp_list);
  infofile << std::endl;

  infofile << "Post-processed data format" << std::endl;
  infofile << "1.  energy" << std::endl;
  infofile << "2.  number of events" << std::endl;
  infofile << "3.  start time [s] and start coordinates [cm]" << std::endl;
  infofile << "4.  energy loss [eV] for each mechanism" << std::endl;
  infofile << "5.  number of ionizations per element" << std::endl;
  infofile << "6.  number of secondary electrons per secondary energy bin" << std::endl;
  infofile << "7.  average parallel distance [cm] for each time bin" << std::endl;
  infofile << "8.  average perpendicular distance [cm] for each time bin" << std::endl;
  infofile << "9.  average energy [eV] for each time bin" << std::endl;
  infofile << "10. energy loss [eV] per time bin" << std::endl;
  infofile << "11. energy loss [eV] per parallel and perpendicular distance bin" << std::endl;

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
void processEvent(const Event* event, std::vector<PartData>& part_data_list, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list) {

  if ( event == nullptr ) {
    std::cerr << "Error: Null event pointer." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

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

    // increment the event count
    part_data.num_ev += 1;
    
    // compute time index
    double t_rel = event->time - part_data.t_start;
    size_t idx_time = findIdx(t_rel, time_list);
    
    // compute distance indices
    double x_rel = event->x - part_data.x_start;
    double y_rel = event->y - part_data.y_start;
    double z_rel = event->z - part_data.z_start;
    double dis_perp = sqrt(x_rel*x_rel + y_rel*y_rel);
    size_t idx_dis_par = findIdx(z_rel, dis_par_list);
    size_t idx_dis_perp = findIdx(dis_perp, dis_perp_list);

    // compute energy loss
    double ener_loss = event->ener_loss + event->ener_loss_moller + event->ener_loss_sync + event->ener_loss_cher;

    // compute time histograms
    if ( idx_time > 0 && idx_time < time_list.size() ) {
      part_data.num_ev_time[idx_time - 1] += 1;
      part_data.dis_par_time[idx_time - 1] += z_rel;
      part_data.dis_perp_time[idx_time - 1] += dis_perp;
      part_data.ener_time[idx_time - 1] += event->ener;
      part_data.ener_loss_time[idx_time - 1] += ener_loss;
    }

    // compute distance histograms
    if ( idx_dis_par > 0 && idx_dis_par < dis_par_list.size() && idx_dis_perp > 0 && idx_dis_perp < dis_perp_list.size() ) {
      part_data.ener_loss_dis2d[idx_dis_par - 1][idx_dis_perp - 1] += ener_loss;
    }
    
    // compute interaction histograms
    size_t idx_ener_sec;
    switch ( event->interaction ) {
      case flags::brem:
      part_data.ener_loss_mech[0] += event->ener_loss;
      break;
      case flags::exc:
      part_data.ener_loss_mech[1] += event->ener_loss;
      break;
      case flags::ion:
      part_data.ener_loss_mech[2] += event->ener_loss;
      part_data.num_ion_elem[event->Zelem - 1] += 1;
      idx_ener_sec = findIdx(event->ener_sec, ener_sec_list);
      if ( idx_ener_sec > 0 && idx_ener_sec < ener_sec_list.size() ) {
        part_data.num_sec_ener[idx_ener_sec - 1] += 1;
      }
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

/**
 * @brief Post-process particle data.
 * 
 * @param part_data_list The list of particle data structs for each energy.
 */
void postProcPartData(std::vector<PartData>& part_data_list) {
  for ( size_t i = 0; i < part_data_list.size(); i++ ) {
    PartData &part_data = part_data_list[i];
    for ( size_t j = 0; j < part_data.num_ev_time.size(); j++ ) {
      if ( part_data.num_ev_time[j] > 0 ) {
        part_data.dis_par_time[j] /= part_data.num_ev_time[j];
        part_data.dis_perp_time[j] /= part_data.num_ev_time[j];
        part_data.ener_time[j] /= part_data.num_ev_time[j];
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
void writePartData(std::string &outfile_name, std::vector<PartData>& part_data_list) {
  
  std::ofstream outfile(outfile_name, std::ios::app);

  if ( !outfile ) {
    std::cerr << "Failed to open file for writing." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  for ( size_t i = 0; i < part_data_list.size(); i++ ) {
    PartData &part_data = part_data_list[i];

    outfile << part_data.ener << std::endl;
    outfile << part_data.num_ev << std::endl;
    outfile << part_data.t_start << ", " << part_data.x_start << "," << part_data.y_start << "," << part_data.z_start << std::endl;
    writeVector(outfile, part_data.ener_loss_mech);
    writeVector(outfile, part_data.num_ion_elem);
    writeVector(outfile, part_data.num_sec_ener);
    writeVector(outfile, part_data.dis_par_time);
    writeVector(outfile, part_data.dis_perp_time);
    writeVector(outfile, part_data.ener_time);
    writeVector(outfile, part_data.ener_loss_time);
    writeVector2d(outfile, part_data.ener_loss_dis2d);
  }

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
void processFile(std::string& datafile_name, std::string& outfile_name, size_t num_event_per_chunk, const std::vector<double> &ener_list, const std::vector<double> &ener_sec_list, const std::vector<double> &time_list, const std::vector<double> &dis_par_list, const std::vector<double> &dis_perp_list) {

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

  size_t npart = 0;
  int current_id = -1;
  std::vector<PartData> part_data_list(ener_list.size(), PartData(ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size()));
  
  while ( datafile.read(buffer.data(), chunk_size) ) {
    for ( size_t i = 0; i < num_event_per_chunk; i++ ) {
      event = reinterpret_cast<Event*>(buffer.data() + i * event_size);
      if ( event->id != current_id ) {
        if ( current_id != -1 ) {
          postProcPartData(part_data_list);
          writePartData(outfile_name, part_data_list);
        }
        for ( size_t j = 0; j < ener_list.size(); j++ ) {
          part_data_list[j].reset(event->id, ener_list[j], ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size());
        }
        current_id = event->id;
        npart++;
      } else {
        processEvent(event, part_data_list, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list);
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
            postProcPartData(part_data_list);
            writePartData(outfile_name, part_data_list);
          }
          for ( size_t j = 0; j < ener_list.size(); j++ ) {
            part_data_list[j].reset(event->id, ener_list[j], ener_sec_list.size(), time_list.size(), dis_par_list.size(), dis_perp_list.size());
          }
          current_id = event->id;
          npart++;
        } else {
          processEvent(event, part_data_list, ener_list, ener_sec_list, time_list, dis_par_list, dis_perp_list);
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
      std::perror("Error deleting file.");
    }
  }

  outfile.close();
  if ( !outfile.good() ) {
    std::cerr << "Error reading file " << outfile_name << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}
