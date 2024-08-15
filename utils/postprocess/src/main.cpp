// includes
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <mpi.h>

// headers
#include "io.h"
#include "functions.h"
#include "parser.h"
#include "const.h"
#include "process.h"

// types
template <typename T>
using vector2d = std::vector<std::vector<T>>;

int main(int argc, char** argv) {

  // initialize MPI
  MPI_Init(&argc, &argv);

  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // start timer
  auto start = std::chrono::steady_clock::now();

  // read config file
  Config config;
  std::string configfile_name = "../config/config.ini";
  parseConfig(configfile_name, config);
  if ( rank == 0 ) std::cout << "Read config file." << std::endl;

  const std::string outfile_name = config["IO"]["outpath"] + "/data.txt";
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);
  const double mach_A = std::stod(config["Bfield"]["mach_A"]);
  const double L = std::stod(config["Bfield"]["L"]);

  // read geometry
  int geo;
  std::string geo_str = std::string(config["Geometry"]["geo"]);
  if ( geo_str == "plane" ) {
    geo = geo_tag::plane;
  } else if ( geo_str == "cylinder" ) {
    geo = geo_tag::cylinder;
  } else {
    geo = geo_tag::none;
  }

  // make lists
  size_t num_ener, num_escape, num_ener_sec, num_time;
  std::vector<double> ener_list, escape_list, ener_sec_list, time_list;
  makeList(config["Bin.Ener"], ener_list, num_ener);
  makeList(config["Bin.Escape"], escape_list, num_escape);
  makeList(config["Bin.EnerSec"], ener_sec_list, num_ener_sec);
  makeList(config["Bin.Time"], time_list, num_time, constants::hr);

  // write info file
  if ( rank == 0 ) {
    const std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, ener_list, escape_list, ener_sec_list);
  }

  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  MPI_Barrier(MPI_COMM_WORLD);
  int count_loc = 0;

  // create a grid of data structs for each energy and escape length
  std::vector<Data> data_list;
  vector2d<Data> data_grid;
  for ( size_t i = 0; i < num_ener; i++ ) {
    data_list.clear();
    for ( size_t j = 0; j < num_escape; j++ ) {
        data_list.push_back(Data(ener_list[i], escape_list[j], mach_A, L, geo, ener_sec_list));
    }
    data_grid.push_back(data_list);
  }

  // process files
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    processFile(datafile_name, num_event_per_chunk, count_loc, data_grid);
  }

  // flatten the data
  std::vector<double> ener_loss_mech_flat_loc, num_ion_elem_flat_loc, num_sec_ener_flat_loc;
  getFlatData(data_grid, ener_loss_mech_flat_loc, num_ion_elem_flat_loc, num_sec_ener_flat_loc);
  size_t ener_loss_mech_flat_size = num_ener * num_escape * num_mech;
  size_t num_ion_elem_size = num_ener * num_escape * num_elem;
  size_t num_sec_ener_flat_size = num_ener * num_escape * num_ener_sec;
  std::vector<double> ener_loss_mech_flat(ener_loss_mech_flat_size);
  std::vector<double> num_ion_elem_flat(num_ion_elem_size);
  std::vector<double> num_sec_ener_flat(num_sec_ener_flat_size);

  // collect the data on rank 0
  int count_glob;
  MPI_Barrier(MPI_COMM_WORLD);
  if ( rank == 0 ) std::cout << "Collecting data on rank 0." << std::endl;
  MPI_Reduce(&count_loc, &count_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(ener_loss_mech_flat_loc.data(), ener_loss_mech_flat.data(), ener_loss_mech_flat_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(num_ion_elem_flat_loc.data(),   num_ion_elem_flat.data(),   num_ion_elem_size,        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(num_sec_ener_flat_loc.data(),   num_sec_ener_flat.data(),   num_sec_ener_flat_size,   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if ( rank == 0 ) {

    // normalize by particle count
    normalize(ener_loss_mech_flat, static_cast<double>(count_glob));
    normalize(num_ion_elem_flat,   static_cast<double>(count_glob));
    normalize(num_sec_ener_flat,   static_cast<double>(count_glob));
    
    // write data
    std::cout << "Writing data to output file." << std::endl << std::endl;
    clearFile(outfile_name);
    writeData(
      outfile_name, 
      ener_list, escape_list, 
      num_ener_sec, 
      ener_loss_mech_flat, num_ion_elem_flat, num_sec_ener_flat
    );
    
    // compute runtime
    auto now = std::chrono::steady_clock::now();
    auto time_pp = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Post processing complete." << std::endl;
    std::cout << "Packet count: " << count_glob << std::endl;
    std::cout << "Runtime [s]:  " << time_pp << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
