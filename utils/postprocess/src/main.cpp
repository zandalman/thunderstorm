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
template <typename T>
using vector3d = std::vector<vector2d<T>>;
template <typename T>
using vector4d = std::vector<vector3d<T>>;

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
  const std::string histdir_name = config["IO"]["outpath"] + "/hist";
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);
  const int num_hist = std::stoi(config["IO"]["num_hist"]);

  // get parameters
  int geo;
  const double L = std::stod(config["Parameters"]["L"]);
  const double ener_min = std::stod(config["Parameters"]["ener_min"]);
  std::string geo_str = std::string(config["Parameters"]["geo"]);
  if ( geo_str == "plane" ) {
    geo = geo_tag::plane;
  } else if ( geo_str == "sphere" ) {
    geo = geo_tag::sphere;
  } else {
    geo = geo_tag::none;
  }

  // make bin list
  size_t num_rho, num_ener, num_escape, num_mach, num_ener_sec, num_time;
  std::vector<double> rho_list, ener_list, escape_list, mach_list, ener_sec_list, time_list;
  makeList(config["Grid.Rho"], rho_list, num_rho);
  makeList(config["Grid.Ener"], ener_list, num_ener);
  makeList(config["Grid.Escape"], escape_list, num_escape);
  makeList(config["Grid.Mach"], mach_list, num_mach);
  makeList(config["Bin.EnerSec"], ener_sec_list, num_ener_sec);
  makeList(config["Bin.Time"], time_list, num_time, constants::hr);
  vector2d<double> bin_list = {rho_list, mach_list, ener_list, escape_list, ener_sec_list, time_list};

  // define statistics
  std::vector<Stat> stat_list;
  stat_list.push_back(Stat(1, "eps_thm", "thermalization efficiency"));
  stat_list.push_back(Stat(num_inter, "num_ev_inter", "number of events for each interaction"));
  stat_list.push_back(Stat(num_mech, "ener_loss_mech", "energy loss [eV] for each mechanism"));
  stat_list.push_back(Stat(num_elem, "num_ion_elem", "number of ionizations per element"));
  stat_list.push_back(Stat(num_ener_sec - 1, "num_sec_ener", "number of secondary electrons per secondary energy bin"));
  stat_list.push_back(Stat(num_time - 1, "ener_loss_time", "energy loss [eV] per time bin"));
  stat_list.push_back(Stat(num_ener - 1, "time_ener", "time [s] spent per energy bin"));

  // write info file
  if ( rank == 0 ) {
    const std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, bin_list, stat_list);
  }

  // create a grid of data structs
  vector4d<Data> data_grid;
  data_grid.resize(num_rho);
  for (size_t i1 = 0; i1 < num_rho; i1++) {
    data_grid[i1].resize(num_mach);
    for (size_t i2 = 0; i2 < num_mach; i2++) {
      data_grid[i1][i2].resize(num_ener);
      for (size_t i3 = 0; i3 < num_ener; i3++) {
        data_grid[i1][i2][i3].resize(num_escape);
        for ( size_t i4 = 0; i4 < num_escape; i4++ ) {
          data_grid[i1][i2][i3][i4] = Data(rho_list[i], mach_list[i], ener_list[j], escape_list[k], ener_min, L, geo, stat_list);
        }
      }
    }
  }

  // get data file indices for this rank
  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  // get number of histories for this rank
  int num_hist_per_rank = static_cast<int>(ceil(num_hist / size));
  int idx_hist_min = std::min(rank * num_hist_per_rank, num_hist);
  int idx_hist_max = std::min((rank + 1) * num_hist_per_rank, num_hist);

  // process files
  MPI_Barrier(MPI_COMM_WORLD);
  int count = 0;
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    processFile(
      datafile_name, 
      num_event_per_chunk, 
      bin_list, 
      stat_list, 
      histdir_name,
      idx_hist_min, 
      idx_hist_max, 
      count, 
      data_grid
    );
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
    calcMoment(size_flat, count, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat, var_stat_list_flat, skew_stat_list_flat, kurt_stat_list_flat);
    
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
