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

  // get file paths
  const std::string outfile_name = config["IO"]["outpath"] + "/data.txt";
  const std::string histdir_name = config["IO"]["outpath"] + "/hist";
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);
  const int num_hist = std::stoi(config["IO"]["num_hist"]);

  // get miscellaneous parameters
  int walltime = std::stoi(config["Misc"]["walltime"]);
  double rho_sim = std::stod(config["Misc"]["rho_sim"]);
  double ener_min = std::stod(config["Misc"]["ener_min"]);
  double vmax = std::stod(config["Misc"]["vmax"]) * constants::c;
  size_t ndim = std::stoul(config["Misc"]["ndim"]);
  size_t nmom = std::stoul(config["Misc"]["nmom"]);

  // make bin list
  size_t num_mach, num_col, num_ener, num_ener_sec;
  std::vector<double> mach_list, col_list, ener_list, ener_sec_list, time_list;
  makeList(config["Grid.Mach"], mach_list, num_mach);
  makeList(config["Grid.Sigma"], col_list, num_col);
  makeList(config["Bin.Ener"], ener_list, num_ener);
  makeList(config["Bin.EnerSec"], ener_sec_list, num_ener_sec);
  vector2d<double> bin_list = {mach_list, col_list, ener_list, ener_sec_list};

  // define statistics
  std::vector<Stat> stat_list;
  stat_list.resize(6);
  stat_list[0] = Stat(1, "ener_thm", "thermalized energy [eV]");
  stat_list[1] = Stat(num_ener - 1, "ener", "energy spectrum of supra-thermal electrons [eV]");
  stat_list[2] = Stat(num_ener_sec - 1, "ener_sec", "energy spectrum of secondary electrons [eV]");
  stat_list[3] = Stat(num_ener - 1, "time_ener", "time [s] spent per energy bin");
  stat_list[4] = Stat(num_inter, "ener_loss_inter", "energy loss by interaction mechanism [eV]");
  stat_list[5] = Stat(num_elem, "num_ion_elem", "number of ionizations per element per ion stage");

  // write info file
  if ( rank == 0 ) {
    const std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, ndim, nmom, bin_list, stat_list);
  }

  // just allocate 3d obj with vol * sizeof(Data)

  // create a grid of data structs
  // vector3d<Data> data_grid(num_mach, vector2d<Data>(num_col, std::vector<Data>(num_ener - 1)));
  vector3d<Data> data_grid;
  data_grid.resize(num_mach);
  for (size_t i = 0; i < num_mach; i++) {
    data_grid[i].resize(num_col);
    for (size_t j = 0; j < num_col; j++) {
      data_grid[i][j].reserve(num_ener - 1);
      for ( size_t k = 0; k < num_ener - 1; k++ ) {
        data_grid[i][j].emplace_back(std::move(Data(
          mach_list[i], 
          col_list[j] / rho_sim, 
          ener_list[k], 
          ener_list[k+1], 
          ener_min,
          col_list[j] / rho_sim / vmax, 
          ndim, 
          nmom,
          stat_list
        )));
      }
    }
  }

  // get data file indices for this rank
  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(static_cast<double>(num_file) / static_cast<double>(size)));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  // get number of histories for this rank
  int num_hist_per_rank = static_cast<int>(ceil(static_cast<double>(num_hist) / static_cast<double>(size)));;
  int idx_hist_min = std::min(rank * num_hist_per_rank, num_hist);
  int idx_hist_max = std::min((rank + 1) * num_hist_per_rank, num_hist);

  // process files
  MPI_Barrier(MPI_COMM_WORLD);
  int count = 0;
  int idx_hist = idx_hist_min;
  bool no_time = false;
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    if ( no_time ) break;
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    processFile(
      datafile_name, 
      num_event_per_chunk,
      start, 
      walltime,
      bin_list, 
      stat_list, 
      histdir_name,
      idx_hist_max, 
      idx_hist,
      count, 
      data_grid,
      no_time
    );
  }

  // flatten the data
  size_t size_flat;
  std::vector<double> M1_stat_list_flat, M2_stat_list_flat, M3_stat_list_flat, M4_stat_list_flat;
  getFlatData(
    data_grid, 
    stat_list, 
    2 * ndim,
    nmom,
    M1_stat_list_flat, 
    M2_stat_list_flat, 
    M3_stat_list_flat, 
    M4_stat_list_flat, 
    size_flat
  );
  
  // collect the data on rank 0
  MPI_Barrier(MPI_COMM_WORLD);
  if ( rank == 0 ) {
    
    std::cout << "Aggregating statistics on rank 0." << std::endl;
    std::cout << "Collected ranks: |";
    
    int count_other;
    std::vector<double> M1_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M2_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M3_stat_list_flat_other(size_flat, 0.);
    std::vector<double> M4_stat_list_flat_other(size_flat, 0.);
    
    for (int i = 1; i < size; i++) {
      MPI_Recv(&count_other, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(M1_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if ( nmom >= 2 ) MPI_Recv(M2_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if ( nmom >= 3 ) MPI_Recv(M3_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if ( nmom >= 4 ) MPI_Recv(M4_stat_list_flat_other.data(), size_flat, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      addStat(
        size_flat, 
        nmom,
        count_other, 
        M1_stat_list_flat_other, 
        M2_stat_list_flat_other, 
        M3_stat_list_flat_other, 
        M4_stat_list_flat_other, 
        count, 
        M1_stat_list_flat, 
        M2_stat_list_flat, 
        M3_stat_list_flat, 
        M4_stat_list_flat
      );
      std::cout << "|";
    } 
    std::cout << std::endl << std::endl;

  } else {
    
    MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(M1_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    if ( nmom >= 2 ) MPI_Send(M2_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    if ( nmom >= 3 ) MPI_Send(M3_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    if ( nmom >= 4 ) MPI_Send(M4_stat_list_flat.data(), size_flat, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
  }

  if ( rank == 0 ) {

    // compute central moments
    std::vector<double> var_stat_list_flat(size_flat, 0.0);
    std::vector<double> skew_stat_list_flat(size_flat, 0.0);
    std::vector<double> kurt_stat_list_flat(size_flat, 0.0);
    if ( nmom >= 2 ) {
      calcMoment(
        size_flat, 
        nmom, 
        count, 
        M2_stat_list_flat, 
        M3_stat_list_flat, 
        M4_stat_list_flat, 
        var_stat_list_flat, 
        skew_stat_list_flat, 
        kurt_stat_list_flat
      );
    }
    
    // write data
    std::cout << "Writing data to output file." << std::endl << std::endl;
    clearFile(outfile_name);
    writeData(
      outfile_name, 
      bin_list, 
      stat_list, 
      2 * ndim,
      nmom,
      M1_stat_list_flat, 
      var_stat_list_flat, 
      skew_stat_list_flat, 
      kurt_stat_list_flat
    );
    
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
