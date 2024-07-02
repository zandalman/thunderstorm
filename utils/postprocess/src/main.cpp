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

  std::cout << rank << ", " << size << std::endl;

  // start timer
  auto start = std::chrono::steady_clock::now();

  // read config file
  std::string config_dir = "../config/config.ini";
  Config config;
  parseConfig(config_dir, config);
  if ( rank == 0 ) std::cout << "Read config file." << std::endl;
  const std::string data_path = config["IO"]["data_path"];
  const size_t num_event_per_chunk = std::stoul(config["IO"]["num_event_per_chunk"]);

  // make energy list
  std::vector<double> ener_list;
  double ener_min = std::stod(config["Energy"]["min"]);
  double ener_max = std::stod(config["Energy"]["max"]);
  size_t num_ener = std::stoul(config["Energy"]["num"]);
  bool log_ener = config["Energy"]["log"] == "true";
  linspace(ener_min, ener_max, num_ener, log_ener, ener_list);

  // make time list
  std::vector<double> time_list;
  double time_min = std::stod(config["Time"]["min"]) * constants::hr;
  double time_max = std::stod(config["Time"]["max"]) * constants::hr;
  size_t num_time = std::stoul(config["Time"]["num"]);
  bool log_time = config["Time"]["log"] == "true";
  linspace(time_min, time_max, num_time, log_time, time_list);

  // make distance list
  std::vector<double> dis_list;
  double dis_min = std::stod(config["Distance"]["min"]) * constants::AU;
  double dis_max = std::stod(config["Distance"]["max"]) * constants::AU;
  size_t num_dis = std::stoul(config["Distance"]["num"]);
  bool log_dis = config["Distance"]["log"] == "true";
  linspace(dis_min, dis_max, num_dis, log_dis, dis_list);

  // write info file
  if ( rank == 0 ) {
    std::string infofile = config["IO"]["outpath"] + "/info.txt";
    writeInfo(infofile, ener_list, time_list, dis_list);
  }

  int num_file = std::stoi(config["IO"]["num_file"]);
  int num_file_per_rank = static_cast<int>(ceil(num_file / size));
  int idx_file_min = std::min(rank * num_file_per_rank, num_file - 1);
  int idx_file_max = std::min((rank + 1) * num_file_per_rank, num_file - 1);

  MPI_Barrier(MPI_COMM_WORLD);
  for ( int i = idx_file_min; i <= idx_file_max; i++ ) {
    std::string filename = data_path + "/data.bin." + std::to_string(i);
    std::string outfile = config["IO"]["outpath"] + "/data.txt." + std::to_string(i);
    clearFile(outfile);
    processFile(filename, outfile, num_event_per_chunk, ener_list, time_list, dis_list);
  }

  if ( rank == 0 ) {
    // concatenateFiles(data_path, num_file);
    auto now = std::chrono::steady_clock::now();
    auto tpp = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
    std::cout << "Post processing complete." << std::endl;
    std::cout << "Runtimes [s]: " << tpp << std::endl;
  }

  MPI_Finalize();
  return 0;
}
