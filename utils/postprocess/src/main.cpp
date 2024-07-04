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
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);

  // make lists
  size_t num_ener, num_ener_sec, num_time, num_dis;
  std::vector<double> ener_list, ener_sec_list, time_list, dis_list;
  makeList(config["Energy"], 1., num_ener, ener_list);
  makeList(config["EnergySec"], 1., num_ener_sec, ener_sec_list);
  makeList(config["Time"], constants::hr, num_time, time_list);
  makeList(config["Distance"], constants::AU, num_dis, dis_list);

  // write info file
  if ( rank == 0 ) {
    std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, config, ener_list, ener_sec_list, time_list, dis_list);
  }

  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file);

  MPI_Barrier(MPI_COMM_WORLD);
  for ( int i = idx_file_min; i < idx_file_max; i++ ) {
    std::string datafile_name = data_path + "/data.bin." + std::to_string(i);
    std::string outfile_name = config["IO"]["outpath"] + "/data.txt." + std::to_string(i);
    clearFile(outfile_name);
    processFile(datafile_name, outfile_name, num_event_per_chunk, ener_list, ener_sec_list, time_list, dis_list);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if ( rank == 0 ) {
    concatenateFiles(config["IO"]["outpath"], num_file);
    auto now = std::chrono::steady_clock::now();
    auto tpp = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Post processing complete." << std::endl;
    std::cout << "Runtime [s]: " << tpp << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
