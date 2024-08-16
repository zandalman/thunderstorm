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

    std::vector<std::string> stat_def_list = {
    "energy loss [eV] for each mechanism",
    "number of ionizations per element",
    "number of secondary electrons per secondary energy bin"
  };

  // read config file
  Config config;
  std::string configfile_name = "../config/config.ini";
  parseConfig(configfile_name, config);
  if ( rank == 0 ) std::cout << "Read config file." << std::endl;

  // get file paths
  const std::string outfile_name = config["IO"]["outpath"] + "/data.txt";
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);
  
  // get magnetic field configuration
  const double mach_A = std::stod(config["Bfield"]["mach_A"]);
  const double L = std::stod(config["Bfield"]["L"]);

  // get geometry
  int geo;
  std::string geo_str = std::string(config["Geometry"]["geo"]);
  if ( geo_str == "plane" ) {
    geo = geo_tag::plane;
  } else if ( geo_str == "cylinder" ) {
    geo = geo_tag::cylinder;
  } else {
    geo = geo_tag::none;
  }

  // make bin list
  size_t num_ener, num_escape, num_ener_sec, num_time;
  std::vector<double> ener_list, escape_list, ener_sec_list, time_list;
  makeList(config["Bin.Ener"], ener_list, num_ener);
  makeList(config["Bin.Escape"], escape_list, num_escape);
  makeList(config["Bin.EnerSec"], ener_sec_list, num_ener_sec);
  makeList(config["Bin.Time"], time_list, num_time, constants::hr);
  vector2d<double> bin_list = {ener_list, escape_list, ener_sec_list, time_list};

  // define statistics
  std::vector<Stat> stat_list;
  stat_list.push_back(Stat(num_mech, "num_ener_mech", "energy loss [eV] for each mechanism"));
  stat_list.push_back(Stat(num_elem, "num_ion_elem", "number of ionizations per element"));
  stat_list.push_back(Stat(num_ener_sec - 1, "num_sec_ener", "number of secondary electrons per secondary energy bin"));
  stat_list.push_back(Stat(num_time - 1, "ener_loss_time", "energy loss [eV] per time bin"));

  // write info file
  if ( rank == 0 ) {
    const std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, ener_list, escape_list, ener_sec_list, stat_list);
  }

  // get data file indices
  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  // create a grid of data structs for each energy and escape length
  std::vector<Data> data_list;
  vector2d<Data> data_grid;
  for ( size_t i = 0; i < num_ener; i++ ) {
    data_list.clear();
    for ( size_t j = 0; j < num_escape; j++ ) {
        data_list.push_back(Data(ener_list[i], escape_list[j], mach_A, L, geo, stat_list));
    }
    data_grid.push_back(data_list);
  }

  // process files
  MPI_Barrier(MPI_COMM_WORLD);
  int count = 0;
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    processFile(datafile_name, num_event_per_chunk, bin_list, stat_list, count, data_grid);
  }

  // flatten the data
  size_t size_flat;
  std::vector<double> mean_stat_list_flat, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat;
  getFlatData(data_grid, stat_list, mean_stat_list_flat, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat, size_flat);
  
  // collect the data on rank 0
  MPI_Barrier(MPI_COMM_WORLD);
  if ( rank == 0 ) {
    
    std::cout << "Aggregating statistics on rank 0." << std::endl;
    std::cout << "Collected ranks: |";
    
    int count_other;
    std::vector<double> mean_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M2_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M3_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M4_stat_list_flat_other(size_flat, 0.);
    
    for (int i = 1; i < size; i++) {
      MPI_Recv(&count_other, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(mean_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(M2_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(M3_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(M4_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      addStat(size_flat, count_other, mean_stat_list_flat_other, M2_stat_list_flat_other, M3_stat_list_flat_other, M4_stat_list_flat_other, count, mean_stat_list_flat, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat);
      std::cout << "|";
    } 
    std::cout << std::endl << std::endl;

  } else {
    
    MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(mean_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(M2_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    MPI_Send(M3_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    MPI_Send(M4_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
  }

  if ( rank == 0 ) {

    // compute central moments
    std::vector<double> var_stat_list_flat(size_flat, 0.0);
    std::vector<double> skew_stat_list_flat(size_flat, 0.0);
    std::vector<double> kurt_stat_list_flat(size_flat, 0.0);
    calcMoment(size_flat, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat, var_stat_list_flat, skew_stat_list_flat, kurt_stat_list_flat);
    
    // write data
    std::cout << "Writing data to output file." << std::endl << std::endl;
    clearFile(outfile_name);
    writeData(outfile_name, bin_list, stat_list, mean_stat_list_flat, var_stat_list_flat, skew_stat_list_flat, kurt_stat_list_flat);
    
    // compute runtime
    auto now = std::chrono::steady_clock::now();
    auto time_pp = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Post processing complete." << std::endl;
    std::cout << "Packet count: " << count << std::endl;
    std::cout << "Runtime [s]:  " << time_pp << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
